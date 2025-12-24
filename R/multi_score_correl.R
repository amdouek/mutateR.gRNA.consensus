#' Summarise agreement between scoring methods (Tiered Analysis)
#'
#' Computes correlation coefficients, rank overlap statistics, and tiered
#' consensus metrics to quantify agreement between on-target scoring methods.
#' Explicitly handles legacy methods (e.g., RuleSet1) that may systematically
#' differ from modern deep learning approaches.
#'
#' @param multi_scores Data frame returned by \code{score_grnas_multi()}.
#' @param methods Character vector. Subset of methods to analyse (default: all).
#' @param modern_methods Character vector. Methods considered "modern" for tiered analysis.
#'        Default: c("deepspcas9", "deephf", "ruleset3") for Cas9,
#'                 c("deepcpf1", "enpamgb") for Cas12a.
#' @param legacy_methods Character vector. Legacy methods to model separately.
#'        Default: c("ruleset1") for Cas9, NULL for Cas12a.
#' @param top_n Integer vector. Thresholds for top-N overlap analysis (default: c(10, 20, 50)).
#' @param cor_method Character. Correlation method: "spearman" (default) or "pearson".
#'
#' @return A list of class "score_agreement" with components:
#'   - correlations: Data frame of pairwise correlation coefficients
#'   - correlation_matrix: Matrix form for heatmap plotting
#'   - rank_overlaps: Data frame of top-N overlap statistics
#'   - summary_stats: Per-method summary statistics
#'   - tiered_consensus: Tiered consensus metrics (modern, robust, legacy divergence
#'   - discordant_examples: gRNAs with largest rank disagreements (modern methods only)
#'   - legacy_divergent: gRNAs where legacy methods most diverge from modern consensus
#'   - concordant_top: gRNAs ranking well across methods
#'   - method_clusters: Hierarchical clustering of methods by correlation
#'   - metadata: Analysis parameters
#'
#' @export
summarize_score_agreement <- function(multi_scores,
                                      methods = NULL,
                                      modern_methods = NULL,
                                      legacy_methods = NULL,
                                      top_n = c(10, 20, 50),
                                      cor_method = c("spearman", "pearson")) {

  cor_method <- match.arg(cor_method)

  # ═══════════════════════════════════════════════════════════════
  # Identify and validate methods
  # ═══════════════════════════════════════════════════════════════

  stored_methods <- attr(multi_scores, "methods")
  nuclease <- attr(multi_scores, "nuclease")

  if (is.null(stored_methods)) {
    meta_cols <- c("grna_id", "protospacer_sequence", "sequence_context",
                   "exon_rank", "chr", "start", "end", "strand")
    stored_methods <- setdiff(names(multi_scores), meta_cols)
    stored_methods <- stored_methods[sapply(multi_scores[stored_methods], is.numeric)]
  }

  if (is.null(methods)) {
    methods <- stored_methods
  }

  # Set default modern/legacy classification based on nuclease
  if (is.null(modern_methods)) {
    if (is.null(nuclease) || nuclease == "Cas9") {
      modern_methods <- intersect(c("deepspcas9", "deephf", "ruleset3"), methods)
    } else {
      modern_methods <- intersect(c("deepcpf1", "enpamgb"), methods)
    }
  }

  if (is.null(legacy_methods)) {
    if (is.null(nuclease) || nuclease == "Cas9") {
      legacy_methods <- intersect(c("ruleset1"), methods)
    } else {
      legacy_methods <- character(0)
    }
  }

  # Validate
  modern_methods <- intersect(modern_methods, methods)
  legacy_methods <- intersect(legacy_methods, methods)

  n_grnas <- nrow(multi_scores)
  score_mat <- multi_scores[, methods, drop = FALSE]

  # ═══════════════════════════════════════════════════════════════
  # 1. Pairwise Correlations
  # ═══════════════════════════════════════════════════════════════

  pairs <- combn(methods, 2, simplify = FALSE)

  cor_results <- lapply(pairs, function(pair) {
    m1 <- pair[1]
    m2 <- pair[2]

    valid <- complete.cases(score_mat[, c(m1, m2)])
    n_valid <- sum(valid)

    if (n_valid < 3) {
      return(data.frame(
        method1 = m1, method2 = m2,
        correlation = NA_real_, p_value = NA_real_, n_valid = n_valid,
        comparison_type = classify_comparison(m1, m2, modern_methods, legacy_methods),
        stringsAsFactors = FALSE
      ))
    }

    test <- cor.test(score_mat[[m1]], score_mat[[m2]],
                     method = cor_method,
                     use = "complete.obs")

    data.frame(
      method1 = m1,
      method2 = m2,
      correlation = test$estimate,
      p_value = test$p.value,
      n_valid = n_valid,
      comparison_type = classify_comparison(m1, m2, modern_methods, legacy_methods),
      stringsAsFactors = FALSE
    )
  })

  correlations_df <- do.call(rbind, cor_results)
  rownames(correlations_df) <- NULL

  # Also create correlation matrix
  cor_matrix <- cor(score_mat, use = "pairwise.complete.obs", method = cor_method)

  # ═══════════════════════════════════════════════════════════════
  # 2. Method Clustering
  # ═══════════════════════════════════════════════════════════════

  method_clusters <- NULL
  if (length(methods) >= 3) {
    tryCatch({
      dist_mat <- as.dist(1 - cor_matrix)
      method_clusters <- hclust(dist_mat, method = "complete")
    }, error = function(e) NULL)
  }

  # ═══════════════════════════════════════════════════════════════
  # 3. Compute Ranks
  # ═══════════════════════════════════════════════════════════════

  ranks <- as.data.frame(lapply(methods, function(m) {
    scores <- score_mat[[m]]
    rank(-scores, na.last = "keep", ties.method = "average")
  }))
  names(ranks) <- methods

  # ═══════════════════════════════════════════════════════════════
  # 4. Tiered Consensus Metrics
  # ═══════════════════════════════════════════════════════════════

  tiered <- data.frame(
    grna_id = multi_scores$grna_id,
    protospacer_sequence = multi_scores$protospacer_sequence,
    stringsAsFactors = FALSE
  )

  # --- Tier 1: Modern Method Consensus ---
  if (length(modern_methods) >= 2) {
    modern_ranks <- ranks[, modern_methods, drop = FALSE]

    tiered$mean_rank_modern <- rowMeans(modern_ranks, na.rm = TRUE)
    tiered$median_rank_modern <- apply(modern_ranks, 1, median, na.rm = TRUE)
    tiered$max_rank_modern <- apply(modern_ranks, 1, max, na.rm = TRUE)
    tiered$min_rank_modern <- apply(modern_ranks, 1, min, na.rm = TRUE)
    tiered$rank_variance_modern <- apply(modern_ranks, 1, var, na.rm = TRUE)
    tiered$rank_range_modern <- tiered$max_rank_modern - tiered$min_rank_modern

    # Quality percentile (modern)
    tiered$quality_percentile_modern <- 100 * (1 - (tiered$mean_rank_modern - 1) / (n_grnas - 1))

    # Concordance score (modern)
    max_var_modern <- max(tiered$rank_variance_modern, na.rm = TRUE)
    tiered$concordance_modern <- if (max_var_modern > 0) {
      100 * (1 - tiered$rank_variance_modern / max_var_modern)
    } else {
      rep(100, n_grnas)
    }
  }

  # --- Tier 2: Robust Full Consensus ---
  tiered$mean_rank_all <- rowMeans(ranks, na.rm = TRUE)
  tiered$median_rank_all <- apply(ranks, 1, median, na.rm = TRUE)
  tiered$trimmed_mean_rank_all <- apply(ranks, 1, mean, trim = 0.25, na.rm = TRUE)
  tiered$max_rank_all <- apply(ranks, 1, max, na.rm = TRUE)
  tiered$rank_variance_all <- apply(ranks, 1, var, na.rm = TRUE)
  tiered$rank_iqr_all <- apply(ranks, 1, IQR, na.rm = TRUE)

  # Quality percentiles
  tiered$quality_percentile_all <- 100 * (1 - (tiered$mean_rank_all - 1) / (n_grnas - 1))
  tiered$quality_percentile_robust <- 100 * (1 - (tiered$median_rank_all - 1) / (n_grnas - 1))

  # --- Tier 3: Legacy Divergence ---
  if (length(legacy_methods) >= 1 && length(modern_methods) >= 2) {
    for (legacy in legacy_methods) {
      legacy_rank <- ranks[[legacy]]
      col_prefix <- gsub("[^a-z0-9]", "_", tolower(legacy))

      # Divergence from modern consensus
      tiered[[paste0(col_prefix, "_divergence")]] <- legacy_rank - tiered$mean_rank_modern
      tiered[[paste0(col_prefix, "_abs_divergence")]] <- abs(legacy_rank - tiered$mean_rank_modern)
      tiered[[paste0(col_prefix, "_divergence_pct")]] <- 100 * (legacy_rank - tiered$mean_rank_modern) / n_grnas

      # Direction: negative means legacy ranks it better (lower) than modern consensus
      tiered[[paste0(col_prefix, "_more_optimistic")]] <- as.integer(legacy_rank < tiered$mean_rank_modern)
    }
  }

  # Add individual ranks
  for (m in methods) {
    tiered[[paste0("rank_", m)]] <- ranks[[m]]
  }

  # ═══════════════════════════════════════════════════════════════
  # 5. Rank Overlap Analysis
  # ═══════════════════════════════════════════════════════════════

  overlap_results <- list()

  for (n in top_n) {
    if (n > n_grnas) next

    for (pair in pairs) {
      m1 <- pair[1]
      m2 <- pair[2]

      top_m1 <- which(ranks[[m1]] <= n)
      top_m2 <- which(ranks[[m2]] <= n)

      overlap <- length(intersect(top_m1, top_m2))
      union_size <- length(union(top_m1, top_m2))
      jaccard <- if (union_size > 0) overlap / union_size else NA_real_

      expected <- n * n / n_grnas
      enrichment <- if (expected > 0) overlap / expected else NA_real_

      overlap_results[[length(overlap_results) + 1]] <- data.frame(
        method1 = m1,
        method2 = m2,
        top_n = n,
        overlap_count = overlap,
        jaccard_index = jaccard,
        expected_overlap = expected,
        enrichment = enrichment,
        comparison_type = classify_comparison(m1, m2, modern_methods, legacy_methods),
        stringsAsFactors = FALSE
      )
    }
  }

  rank_overlaps_df <- do.call(rbind, overlap_results)
  rownames(rank_overlaps_df) <- NULL

  # ═══════════════════════════════════════════════════════════════
  # 6. Per-Method Summary Statistics
  # ═══════════════════════════════════════════════════════════════

  summary_stats <- lapply(methods, function(m) {
    scores <- score_mat[[m]]
    valid_scores <- scores[!is.na(scores)]

    data.frame(
      method = m,
      method_type = if (m %in% modern_methods) "modern" else if (m %in% legacy_methods) "legacy" else "other",
      n_scored = length(valid_scores),
      n_na = sum(is.na(scores)),
      min = min(valid_scores, na.rm = TRUE),
      q25 = quantile(valid_scores, 0.25, na.rm = TRUE),
      median = median(valid_scores, na.rm = TRUE),
      mean = mean(valid_scores, na.rm = TRUE),
      q75 = quantile(valid_scores, 0.75, na.rm = TRUE),
      max = max(valid_scores, na.rm = TRUE),
      sd = sd(valid_scores, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary_stats_df <- do.call(rbind, summary_stats)
  rownames(summary_stats_df) <- NULL

  # ═══════════════════════════════════════════════════════════════
  # 7. Identify Discordant Examples (Modern Methods Only)
  # ═══════════════════════════════════════════════════════════════

  discordant_examples <- NULL

  if (length(modern_methods) >= 2 && "rank_variance_modern" %in% names(tiered)) {
    discordant_df <- tiered[, c("grna_id", "protospacer_sequence",
                                paste0("rank_", modern_methods),
                                "mean_rank_modern", "rank_variance_modern")]
    discordant_df <- discordant_df[order(-discordant_df$rank_variance_modern), ]
    discordant_examples <- head(discordant_df, 20)
  }

  # ═══════════════════════════════════════════════════════════════
  # 8. Legacy Divergent Examples
  # ═══════════════════════════════════════════════════════════════

  legacy_divergent <- NULL

  if (length(legacy_methods) >= 1 && length(modern_methods) >= 2) {
    legacy <- legacy_methods[1]
    col_prefix <- gsub("[^a-z0-9]", "_", tolower(legacy))
    abs_div_col <- paste0(col_prefix, "_abs_divergence")
    div_col <- paste0(col_prefix, "_divergence")

    if (abs_div_col %in% names(tiered)) {
      legacy_div_df <- tiered[, c("grna_id", "protospacer_sequence",
                                  paste0("rank_", legacy),
                                  "mean_rank_modern",
                                  div_col, abs_div_col)]

      # Split into legacy-optimistic (legacy ranks higher) and legacy-pessimistic
      legacy_div_df$direction <- ifelse(tiered[[div_col]] < 0,
                                        "legacy_optimistic",
                                        "legacy_pessimistic")

      legacy_div_df <- legacy_div_df[order(-legacy_div_df[[abs_div_col]]), ]
      legacy_divergent <- head(legacy_div_df, 20)
    }
  }

  # ═══════════════════════════════════════════════════════════════
  # 9. Concordant High Performers
  # ═══════════════════════════════════════════════════════════════

  # Use modern methods if available, otherwise all
  if (length(modern_methods) >= 2 && "max_rank_modern" %in% names(tiered)) {
    concordant_df <- tiered[, c("grna_id", "protospacer_sequence",
                                paste0("rank_", modern_methods),
                                "mean_rank_modern", "max_rank_modern")]
    concordant_df <- concordant_df[order(concordant_df$max_rank_modern), ]
  } else {
    concordant_df <- tiered[, c("grna_id", "protospacer_sequence",
                                paste0("rank_", methods),
                                "mean_rank_all", "max_rank_all")]
    concordant_df <- concordant_df[order(concordant_df$max_rank_all), ]
  }

  concordant_top <- head(concordant_df, 20)

  # ═══════════════════════════════════════════════════════════════
  # 10. Summary Statistics for Tiers
  # ═══════════════════════════════════════════════════════════════

  tier_summary <- list()

  # Modern method correlations
  if (length(modern_methods) >= 2) {
    modern_cors <- correlations_df$correlation[correlations_df$comparison_type == "modern-modern"]
    tier_summary$modern <- list(
      methods = modern_methods,
      mean_correlation = mean(modern_cors, na.rm = TRUE),
      min_correlation = min(modern_cors, na.rm = TRUE),
      max_correlation = max(modern_cors, na.rm = TRUE)
    )
  }

  # Legacy-modern correlations
  if (length(legacy_methods) >= 1 && length(modern_methods) >= 1) {
    legacy_cors <- correlations_df$correlation[correlations_df$comparison_type == "legacy-modern"]
    tier_summary$legacy_modern <- list(
      mean_correlation = mean(legacy_cors, na.rm = TRUE),
      min_correlation = min(legacy_cors, na.rm = TRUE),
      max_correlation = max(legacy_cors, na.rm = TRUE)
    )
  }

  # ═══════════════════════════════════════════════════════════════
  # Return Results
  # ═══════════════════════════════════════════════════════════════

  result <- list(
    correlations = correlations_df,
    correlation_matrix = cor_matrix,
    rank_overlaps = rank_overlaps_df,
    summary_stats = summary_stats_df,
    tiered_consensus = tiered,
    tier_summary = tier_summary,
    discordant_examples = discordant_examples,
    legacy_divergent = legacy_divergent,
    concordant_top = concordant_top,
    method_clusters = method_clusters,
    metadata = list(
      n_grnas = n_grnas,
      methods = methods,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      nuclease = nuclease,
      cor_method = cor_method,
      top_n_thresholds = top_n
    )
  )

  class(result) <- c("score_agreement", class(result))
  return(result)
}


#' Classify comparison type between two methods
#' @noRd
classify_comparison <- function(m1, m2, modern_methods, legacy_methods) {
  m1_type <- if (m1 %in% modern_methods) "modern" else if (m1 %in% legacy_methods) "legacy" else "other"
  m2_type <- if (m2 %in% modern_methods) "modern" else if (m2 %in% legacy_methods) "legacy" else "other"

  types <- sort(c(m1_type, m2_type))
  paste(types, collapse = "-")
}


#' Extract tiered consensus metrics from agreement summary
#'
#' Convenience function to extract the tiered consensus data frame
#' with selected columns for downstream analysis.
#'
#' @param agreement Output from summarize_score_agreement()
#' @param tier Character. One of "modern", "robust", "all", or "legacy"
#' @return Data frame with relevant metrics
#' @export
extract_tiered_metrics <- function(agreement, tier = c("modern", "robust", "all", "legacy")) {
  tier <- match.arg(tier)
  tc <- agreement$tiered_consensus

  base_cols <- c("grna_id", "protospacer_sequence")

  tier_cols <- switch(tier,
                      "modern" = c("mean_rank_modern", "max_rank_modern", "rank_variance_modern",
                                   "quality_percentile_modern", "concordance_modern"),
                      "robust" = c("median_rank_all", "trimmed_mean_rank_all", "rank_iqr_all",
                                   "quality_percentile_robust"),
                      "all" = c("mean_rank_all", "max_rank_all", "rank_variance_all",
                                "quality_percentile_all"),
                      "legacy" = grep("_divergence|_optimistic", names(tc), value = TRUE)
  )

  available_cols <- intersect(c(base_cols, tier_cols), names(tc))
  tc[, available_cols, drop = FALSE]
}
