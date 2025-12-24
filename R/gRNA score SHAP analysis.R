#' Extract sequence features for SHAP analysis
#'
#' Computes a comprehensive set of sequence-based features for each gRNA
#' to enable Shapley value analysis of scoring method behavior.
#'
#' @param scores_df Data frame from score_grnas_multi()
#' @return Data frame with original data plus feature columns
#' @export
extract_grna_features <- function(scores_df) {

  seqs <- toupper(scores_df$protospacer_sequence)
  ctx <- toupper(scores_df$sequence_context)
  n <- length(seqs)

  features <- data.frame(
    grna_id = scores_df$grna_id,
    stringsAsFactors = FALSE
  )

  # --- 1. Global GC composition ---
  features$gc_content <- sapply(seqs, function(s) {
    (str_count(s, "G") + str_count(s, "C")) / nchar(s)
  })

  features$gc_count <- sapply(seqs, function(s) {
    str_count(s, "G") + str_count(s, "C")
  })

  # --- 2. Positional GC (PAM-proximal vs distal) ---
  features$gc_proximal <- sapply(seqs, function(s) {
    # Last 10 nt (PAM-proximal for Cas9)
    prox <- substr(s, nchar(s) - 9, nchar(s))
    (str_count(prox, "G") + str_count(prox, "C")) / 10
  })

  features$gc_distal <- sapply(seqs, function(s) {
    # First 10 nt (PAM-distal)
    dist <- substr(s, 1, 10)
    (str_count(dist, "G") + str_count(dist, "C")) / 10
  })

  # --- 3. Homopolymer runs (capped at 4 to accommodate filtering of runs >4 in mutateR) ---
  features$max_homopolymer <- sapply(seqs, function(s) {
    runs <- gregexpr("(.)\\1+", s, perl = TRUE)[[1]]
    if (runs[1] == -1) return(1L)
    max(attr(runs, "match.length"))
  })

  features$has_polyT <- as.integer(grepl("TTTT", seqs))
  features$has_polyA <- as.integer(grepl("AAAA", seqs))
  features$has_polyG <- as.integer(grepl("GGGG", seqs))
  features$has_polyC <- as.integer(grepl("CCCC", seqs))

  # --- 4. Dinucleotide frequencies ---
  dinucs <- c("AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC",
              "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC")

  for (di in dinucs) {
    features[[paste0("freq_", di)]] <- sapply(seqs, function(s) {
      str_count(s, di) / (nchar(s) - 1)
    })
  }

  # --- 5. Position-specific nucleotides ---
  # Positions known to be important are 1, 2, 3, 17, 18, 19, 20 (Cas9), but for completeness featurise every position.
  key_positions <- c(1:20)

  for (pos in key_positions) {
    for (nuc in c("A", "T", "G", "C")) {
      features[[paste0("pos", pos, "_", nuc)]] <- sapply(seqs, function(s) {
        as.integer(substr(s, pos, pos) == nuc)
      })
    }
  }

  # --- 6. Thermodynamic proxies (placeholder - to populate with additional metrics) ---
  # Simple approximation without external tools
  features$estimated_tm <- sapply(seqs, function(s) {
    gc <- str_count(s, "[GC]")
    at <- str_count(s, "[AT]")
    # Wallace rule approximation (TO DO -- implement Santa Lucia)
    2 * at + 4 * gc
  })

  # --- 7. Self-complementarity (simple) ---
  features$self_complement_score <- sapply(seqs, function(s) {
    # Count potential self-binding dinucleotides
    rc <- chartr("ATGC", "TACG", s)
    rc <- paste(rev(strsplit(rc, "")[[1]]), collapse = "")
    sum(sapply(1:(nchar(s)-1), function(i) {
      substr(s, i, i+1) == substr(rc, i, i+1)
    }))
  })

  # --- 8. Sequence complexity ---
  features$linguistic_complexity <- sapply(seqs, function(s) {
    # Number of unique k-mers vs expected
    kmers_2 <- sapply(1:(nchar(s)-1), function(i) substr(s, i, i+1))
    kmers_3 <- sapply(1:(nchar(s)-2), function(i) substr(s, i, i+2))
    (length(unique(kmers_2)) + length(unique(kmers_3))) /
      (min(16, nchar(s)-1) + min(64, nchar(s)-2))
  })

  # --- 9. PAM-adjacent context ---
  # For Cas9: nucleotide at position -1 (just before PAM)
  features$pre_pam_nuc <- substr(seqs, nchar(seqs), nchar(seqs))
  for (nuc in c("A", "T", "G", "C")) {
    features[[paste0("pre_pam_", nuc)]] <- as.integer(features$pre_pam_nuc == nuc)
  }
  features$pre_pam_nuc <- NULL  # Remove string column

  # --- 10. Exon position ---
  features$exon_rank <- scores_df$exon_rank

  return(features)
}


#' Compute target variables for SHAP analysis (Tiered Version)
#'
#' Calculates outcome variables with explicit handling of modern vs. legacy
#' methods. Provides tiered consensus metrics for flexible analysis that
#' prevents legacy methods from dominating variance-based statistics.
#'
#' @param scores_df Data frame from score_grnas_multi()
#' @param modern_methods Character vector. Modern scoring methods.
#'        Default: c("deepspcas9", "deephf", "ruleset3")
#' @param legacy_methods Character vector. Legacy methods to model separately.
#'        Default: c("ruleset1")
#' @return Data frame with tiered target variables
#' @export
compute_shap_targets <- function(scores_df,
                                 modern_methods = NULL,
                                 legacy_methods = NULL) {

  # Get method metadata
  all_available <- attr(scores_df, "methods")
  nuclease <- attr(scores_df, "nuclease")

  if (is.null(all_available)) {
    meta_cols <- c("grna_id", "protospacer_sequence", "sequence_context",
                   "exon_rank", "chr", "start", "end", "strand")
    all_available <- setdiff(names(scores_df), meta_cols)
    all_available <- all_available[sapply(scores_df[all_available], is.numeric)]
  }

  # Set defaults based on nuclease
  if (is.null(modern_methods)) {
    if (is.null(nuclease) || nuclease == "Cas9") {
      modern_methods <- intersect(c("deepspcas9", "deephf", "ruleset3"), all_available)
    } else {
      modern_methods <- intersect(c("deepcpf1", "enpamgb"), all_available)
    }
  }

  if (is.null(legacy_methods)) {
    if (is.null(nuclease) || nuclease == "Cas9") {
      legacy_methods <- intersect(c("ruleset1"), all_available)
    } else {
      legacy_methods <- character(0)
    }
  }

  all_methods <- union(modern_methods, legacy_methods)
  n_grnas <- nrow(scores_df)

  targets <- data.frame(
    grna_id = scores_df$grna_id,
    stringsAsFactors = FALSE
  )

  # ═══════════════════════════════════════════════════════════════
  # TIER 1: Modern Method Consensus
  # ═══════════════════════════════════════════════════════════════

  if (length(modern_methods) >= 2) {
    modern_scores <- scores_df[, modern_methods, drop = FALSE]

    # Normalise (0-1)
    modern_norm <- as.data.frame(lapply(modern_scores, function(x) {
      rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      if (rng == 0) return(rep(0.5, length(x)))
      (x - min(x, na.rm = TRUE)) / rng
    }))

    # Ranks (1 = best)
    modern_ranks <- as.data.frame(lapply(modern_scores, function(x) {
      rank(-x, na.last = "keep", ties.method = "average")
    }))

    targets$mean_score_modern <- rowMeans(modern_norm, na.rm = TRUE)
    targets$mean_rank_modern <- rowMeans(modern_ranks, na.rm = TRUE)
    targets$max_rank_modern <- apply(modern_ranks, 1, max, na.rm = TRUE)
    targets$min_rank_modern <- apply(modern_ranks, 1, min, na.rm = TRUE)
    targets$rank_variance_modern <- apply(modern_ranks, 1, var, na.rm = TRUE)
    targets$rank_range_modern <- targets$max_rank_modern - targets$min_rank_modern

    # Quality percentile (higher = better)
    targets$quality_percentile_modern <- 100 * (1 - (targets$mean_rank_modern - 1) / (n_grnas - 1))

    # Inverted rank for intuitive interpretation
    targets$quality_inverted_rank_modern <- (n_grnas + 1) - targets$mean_rank_modern

    # Concordance (higher = better agreement)
    max_var <- max(targets$rank_variance_modern, na.rm = TRUE)
    targets$concordance_modern <- if (max_var > 0) {
      100 * (1 - targets$rank_variance_modern / max_var)
    } else {
      rep(100, n_grnas)
    }

    # Discordance (higher = more disagreement)
    targets$discordance_modern <- if (max_var > 0) {
      100 * (targets$rank_variance_modern / max_var)
    } else {
      rep(0, n_grnas)
    }
  }

  # ═══════════════════════════════════════════════════════════════
  # TIER 2: Robust Full Consensus (All Methods)
  # ═══════════════════════════════════════════════════════════════

  all_scores <- scores_df[, all_methods, drop = FALSE]

  all_norm <- as.data.frame(lapply(all_scores, function(x) {
    rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / rng
  }))

  all_ranks <- as.data.frame(lapply(all_scores, function(x) {
    rank(-x, na.last = "keep", ties.method = "average")
  }))

  # Robust statistics (resistant to single outlier)
  targets$median_rank_all <- apply(all_ranks, 1, median, na.rm = TRUE)
  targets$trimmed_mean_rank_all <- apply(all_ranks, 1, mean, trim = 0.25, na.rm = TRUE)
  targets$rank_iqr_all <- apply(all_ranks, 1, IQR, na.rm = TRUE)
  targets$max_rank_all <- apply(all_ranks, 1, max, na.rm = TRUE)

  # Non-robust (for comparison)
  targets$mean_rank_all <- rowMeans(all_ranks, na.rm = TRUE)
  targets$rank_variance_all <- apply(all_ranks, 1, var, na.rm = TRUE)

  # Quality percentiles
  targets$quality_percentile_all <- 100 * (1 - (targets$mean_rank_all - 1) / (n_grnas - 1))
  targets$quality_percentile_robust <- 100 * (1 - (targets$median_rank_all - 1) / (n_grnas - 1))

  # ═══════════════════════════════════════════════════════════════
  # TIER 3: Legacy Divergence
  # ═══════════════════════════════════════════════════════════════

  if (length(legacy_methods) >= 1 && length(modern_methods) >= 2) {
    for (legacy in legacy_methods) {
      legacy_rank <- all_ranks[[legacy]]
      col_prefix <- gsub("[^a-z0-9]", "_", tolower(legacy))

      # Divergence: positive = legacy ranks it worse than modern consensus
      targets[[paste0(col_prefix, "_divergence")]] <- legacy_rank - targets$mean_rank_modern

      # Absolute divergence (magnitude regardless of direction)
      targets[[paste0(col_prefix, "_abs_divergence")]] <- abs(targets[[paste0(col_prefix, "_divergence")]])

      # Scaled divergence (percentage of total gRNAs)
      targets[[paste0(col_prefix, "_divergence_pct")]] <- 100 * targets[[paste0(col_prefix, "_divergence")]] / n_grnas

      # Binary: is legacy more optimistic?
      targets[[paste0(col_prefix, "_more_optimistic")]] <- as.integer(legacy_rank < targets$mean_rank_modern)
    }
  }

  # ═══════════════════════════════════════════════════════════════
  # TIER 4: Binary Selection Flags
  # ═══════════════════════════════════════════════════════════════

  if ("max_rank_modern" %in% names(targets)) {
    targets$is_top20_modern <- as.integer(targets$max_rank_modern <= 20)
    targets$is_top50_modern <- as.integer(targets$max_rank_modern <= 50)
    targets$is_top100_modern <- as.integer(targets$max_rank_modern <= 100)
  }

  targets$is_top50_robust <- as.integer(targets$median_rank_all <= 50)

  # High discordance flags
  if ("rank_variance_modern" %in% names(targets)) {
    targets$is_discordant_modern <- as.integer(
      targets$rank_variance_modern > quantile(targets$rank_variance_modern, 0.9, na.rm = TRUE)
    )
  }

  if (length(legacy_methods) >= 1) {
    legacy <- legacy_methods[1]
    col_prefix <- gsub("[^a-z0-9]", "_", tolower(legacy))
    abs_div_col <- paste0(col_prefix, "_abs_divergence")

    if (abs_div_col %in% names(targets)) {
      targets$is_legacy_divergent <- as.integer(
        targets[[abs_div_col]] > quantile(targets[[abs_div_col]], 0.9, na.rm = TRUE)
      )
    }
  }

  # ═══════════════════════════════════════════════════════════════
  # Individual Method Scores/Ranks
  # ═══════════════════════════════════════════════════════════════

  for (m in all_methods) {
    targets[[paste0("score_", m)]] <- all_scores[[m]]
    targets[[paste0("rank_", m)]] <- all_ranks[[m]]
  }

  # Add metadata as attribute
  attr(targets, "modern_methods") <- modern_methods
  attr(targets, "legacy_methods") <- legacy_methods
  attr(targets, "n_grnas") <- n_grnas

  return(targets)
}

#' Run tiered SHAP analysis across consensus metrics
#'
#' Executes SHAP analysis for multiple target variables representing different
#' aspects of gRNA quality and method agreement. This enables systematic
#' comparison of feature importance across:
#' - Modern method consensus (primary quality metric)
#' - Modern method disagreement (what causes DeepSpCas9/DeepHF/RuleSet3 to differ)
#' - Legacy divergence (what RuleSet1 misses that modern methods catch)
#' - Robust overall quality (outlier-resistant consensus)
#'
#' @param batch_result A mutateR_consensus_batch object
#' @param modern_methods Character vector. Modern scoring methods.
#' @param legacy_methods Character vector. Legacy methods.
#' @param n_top_features Integer. Number of top features to extract (default 25)
#' @param include_method_specific Logical. Also run SHAP for individual methods (default FALSE)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return A list of class "mutateR_tiered_shap" containing:
#'   - modern_quality: SHAP for quality_percentile_modern
#'   - modern_discord: SHAP for discordance_modern
#'   - legacy_divergence: SHAP for legacy divergence
#'   - robust_quality: SHAP for quality_percentile_robust
#'   - comparison: Feature ranking comparison across analyses
#'   - method_specific: (optional) Individual method SHAP results
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Run batch analysis
#' batch <- analyze_score_consensus(
#'   c("TP53", "BRCA1", "EGFR"),
#'   "hsapiens",
#'   BSgenome.Hsapiens.UCSC.hg38
#' )
#'
#' # Run tiered SHAP
#' shap_results <- run_tiered_shap_analysis(batch)
#'
#' # View results
#' print(shap_results)
#'
#' # Access individual analyses
#' shap_results$modern_quality$plots$beeswarm
#' shap_results$legacy_divergence$plots$beeswarm
#'
#' # Compare feature rankings
#' shap_results$comparison
#' }
#'
#' @export
run_tiered_shap_analysis <- function(batch_result,
                                     modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                     legacy_methods = c("ruleset1"),
                                     n_top_features = 25,
                                     include_method_specific = FALSE,
                                     quiet = FALSE) {

  if (!inherits(batch_result, "mutateR_consensus_batch") &&
      !inherits(batch_result, "mutateR_consensus_analysis")) {
    stop("Input must be a mutateR_consensus_batch or mutateR_consensus_analysis object")
  }

  # Convert single analysis to batch format
  if (inherits(batch_result, "mutateR_consensus_analysis")) {
    gene_name <- batch_result$metadata$gene_symbol %||% batch_result$metadata$gene_id
    batch_result <- structure(
      list(batch_result),
      names = gene_name,
      class = c("mutateR_consensus_batch", "list"),
      metadata = list(
        gene_ids = gene_name,
        n_successful = 1,
        methods = batch_result$metadata$methods
      )
    )
  }

  if (!quiet) {
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════════╗\n")
    cat("║           mutateR Tiered SHAP Analysis                        ║\n")
    cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  }

  results <- list()

  # Define target analyses
  tier_targets <- list(
    modern_quality = list(
      target = "quality_percentile_modern",
      description = "Modern method consensus quality",
      interpretation = "Higher = better predicted efficacy by modern methods"
    ),
    modern_discord = list(
      target = "discordance_modern",
      description = "Modern method disagreement",
      interpretation = "Higher = more disagreement among DeepSpCas9/DeepHF/RuleSet3"
    ),
    legacy_divergence = list(
      target = "ruleset1_divergence",
      description = "Legacy (RuleSet1) divergence from modern consensus",
      interpretation = "Positive = RuleSet1 ranks it worse than modern methods"
    ),
    robust_quality = list(
      target = "quality_percentile_robust",
      description = "Robust overall quality (median-based)",
      interpretation = "Higher = better across all methods (outlier-resistant)"
    )
  )

  # Run each tier
  for (tier_name in names(tier_targets)) {
    tier_info <- tier_targets[[tier_name]]

    if (!quiet) {
      cat(sprintf("[%s] %s\n", tier_name, tier_info$description))
      cat(sprintf("    Target: %s\n", tier_info$target))
    }

    results[[tier_name]] <- tryCatch({
      res <- analyze_feature_importance(
        batch_result,
        target = tier_info$target,
        modern_methods = modern_methods,
        legacy_methods = legacy_methods,
        n_top_features = n_top_features
      )
      res$interpretation <- tier_info$interpretation
      res$tier_name <- tier_name

      if (!quiet) cat("    ✓ Complete\n\n")
      res
    }, error = function(e) {
      if (!quiet) cat(sprintf("    ✗ Failed: %s\n\n", e$message))
      NULL
    })
  }

  # Optional: Method-specific analyses
  if (include_method_specific) {
    if (!quiet) cat("Running method-specific analyses...\n")

    results$method_specific <- list()
    all_methods <- union(modern_methods, legacy_methods)

    for (method in all_methods) {
      target_name <- paste0("score_", method)

      if (!quiet) cat(sprintf("  [%s]\n", method))

      results$method_specific[[method]] <- tryCatch({
        analyze_feature_importance(
          batch_result,
          target = target_name,
          modern_methods = modern_methods,
          legacy_methods = legacy_methods,
          n_top_features = n_top_features
        )
      }, error = function(e) {
        if (!quiet) cat(sprintf("    Failed: %s\n", e$message))
        NULL
      })
    }
  }

  # Compare feature rankings across tiers
  if (!quiet) cat("Comparing feature rankings across tiers...\n")
  results$comparison <- compare_tier_rankings(results, n_top = n_top_features)

  # Store metadata
  results$metadata <- list(
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    tiers_completed = names(Filter(Negate(is.null), results[names(tier_targets)])),
    timestamp = Sys.time()
  )

  class(results) <- c("mutateR_tiered_shap", class(results))

  if (!quiet) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("Tiered SHAP analysis complete.\n")
    cat("Access results via: $modern_quality, $modern_discord,\n")
    cat("                    $legacy_divergence, $robust_quality\n")
    cat("Compare features via: $comparison\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")
  }

  return(results)
}


#' Compare feature rankings across SHAP tiers
#' @noRd
compare_tier_rankings <- function(shap_results, n_top = 25) {

  tier_names <- c("modern_quality", "modern_discord", "legacy_divergence", "robust_quality")
  available_tiers <- intersect(tier_names, names(shap_results))
  available_tiers <- available_tiers[sapply(available_tiers, function(t) !is.null(shap_results[[t]]))]

  if (length(available_tiers) < 2) {
    return(NULL)
  }

  # Extract feature importance from each tier
  importance_list <- lapply(available_tiers, function(tier) {
    fi <- shap_results[[tier]]$feature_importance
    if (is.null(fi)) return(NULL)

    # Handle different formats from shapviz
    if (is.data.frame(fi)) {
      if ("feature" %in% names(fi)) {
        features <- fi$feature
      } else {
        features <- rownames(fi)
      }

      # Get importance values
      if ("importance" %in% names(fi)) {
        values <- fi$importance
      } else if ("mean_abs_shap" %in% names(fi)) {
        values <- fi$mean_abs_shap
      } else {
        values <- fi[[1]]
      }

      data.frame(
        feature = features[1:min(n_top, length(features))],
        rank = 1:min(n_top, length(features)),
        tier = tier,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })

  importance_list <- Filter(Negate(is.null), importance_list)

  if (length(importance_list) == 0) return(NULL)

  # Combine all rankings
  all_rankings <- do.call(rbind, importance_list)

  # Get all unique features
  all_features <- unique(all_rankings$feature)

  # Create wide-format comparison
  comparison <- data.frame(feature = all_features, stringsAsFactors = FALSE)

  for (tier in available_tiers) {
    tier_data <- all_rankings[all_rankings$tier == tier, ]
    comparison[[paste0("rank_", tier)]] <- match(comparison$feature, tier_data$feature)
  }

  # Calculate summary statistics
  rank_cols <- grep("^rank_", names(comparison), value = TRUE)

  comparison$mean_rank <- rowMeans(comparison[, rank_cols], na.rm = TRUE)
  comparison$n_top_tiers <- rowSums(!is.na(comparison[, rank_cols]))
  comparison$rank_variance <- apply(comparison[, rank_cols], 1, var, na.rm = TRUE)

  # Sort by mean rank (features important across multiple tiers first)
  comparison <- comparison[order(comparison$mean_rank), ]

  # Add interpretation column
  comparison$interpretation <- sapply(1:nrow(comparison), function(i) {
    row <- comparison[i, ]

    if (row$n_top_tiers == length(available_tiers)) {
      return("UNIVERSAL: Important across all analyses")
    }

    # Check which tiers it's important in
    important_in <- available_tiers[!is.na(unlist(row[rank_cols]))]

    if (length(important_in) == 1) {
      return(paste0("TIER-SPECIFIC: ", important_in))
    }

    if (all(c("modern_quality", "robust_quality") %in% important_in) &&
        !("legacy_divergence" %in% important_in)) {
      return("QUALITY DRIVER: Affects quality but not legacy divergence")
    }

    if ("legacy_divergence" %in% important_in &&
        !("modern_quality" %in% important_in)) {
      return("LEGACY BLIND SPOT: RuleSet1 misses this")
    }

    return(paste0("MIXED: ", paste(important_in, collapse = ", ")))
  })

  rownames(comparison) <- NULL
  return(comparison)
}


#' Print method for tiered SHAP analysis
#' @export
print.mutateR_tiered_shap <- function(x, ...) {

  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║           mutateR Tiered SHAP Analysis Results                ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

  meta <- x$metadata

  cat("Modern methods:  ", paste(meta$modern_methods, collapse = ", "), "\n")
  cat("Legacy methods:  ", paste(meta$legacy_methods, collapse = ", "), "\n")
  cat("Tiers completed: ", paste(meta$tiers_completed, collapse = ", "), "\n\n")

  # Print summary for each tier
  tier_info <- list(
    modern_quality = "Quality (Modern Consensus)",
    modern_discord = "Disagreement (Modern Methods)",
    legacy_divergence = "Legacy Divergence",
    robust_quality = "Quality (Robust/All Methods)"
  )

  for (tier in meta$tiers_completed) {
    if (!is.null(x[[tier]])) {
      cat(sprintf("─── %s ───\n", tier_info[[tier]]))
      cat(sprintf("    gRNAs: %d | Features: %d\n",
                  x[[tier]]$metadata$n_gRNAs,
                  x[[tier]]$metadata$n_features))

      # Top 3 features
      fi <- x[[tier]]$feature_importance
      if (!is.null(fi) && nrow(fi) >= 3) {
        cat("    Top 3 features:\n")
        for (i in 1:3) {
          feature_name <- if ("feature" %in% names(fi)) fi$feature[i] else rownames(fi)[i]
          cat(sprintf("      %d. %s\n", i, feature_name))
        }
      }
      cat("\n")
    }
  }

  # Print comparison highlights
  if (!is.null(x$comparison)) {
    cat("─── Cross-Tier Feature Comparison ───\n")

    # Universal features
    universal <- x$comparison$feature[grepl("UNIVERSAL", x$comparison$interpretation)]
    if (length(universal) > 0) {
      cat("Universal features (important across all tiers):\n")
      cat("  ", paste(head(universal, 5), collapse = ", "), "\n")
    }

    # Legacy blind spots
    blind_spots <- x$comparison$feature[grepl("LEGACY BLIND SPOT", x$comparison$interpretation)]
    if (length(blind_spots) > 0) {
      cat("Legacy blind spots (RuleSet1 misses these):\n")
      cat("  ", paste(head(blind_spots, 5), collapse = ", "), "\n")
    }

    cat("\n")
  }

  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Access individual tier results: $modern_quality, $modern_discord, etc.\n")
  cat("Access feature comparison: $comparison\n")
  cat("Access plots: $<tier>$plots$beeswarm or $<tier>$plots$bar\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  invisible(x)
}

#' Run SHAP analysis for gRNA features
#'
#' Uses a gradient boosting model with SHAP values to identify which sequence
#' features drive high scores and cross-method agreement/disagreement.
#' Supports tiered target variables for nuanced analysis.
#'
#' @param batch_result A mutateR_consensus_batch or mutateR_consensus_analysis object
#' @param target Character. Target variable to model. Common options:
#'        - "quality_percentile_modern": Modern method consensus (recommended primary)
#'        - "quality_percentile_robust": Outlier-resistant quality
#'        - "discordance_modern": Modern method disagreement
#'        - "ruleset1_divergence": Legacy method divergence
#'        - "score_<method>": Individual method score
#' @param modern_methods Character vector of modern scoring methods
#' @param legacy_methods Character vector of legacy scoring methods
#' @param n_top_features Integer. Number of top features to display (default 20)
#' @param xgb_params List. XGBoost parameters (optional, uses sensible defaults)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return List of class "mutateR_shap_analysis" containing model, SHAP values, and plots
#' @export
analyze_feature_importance <- function(batch_result,
                                       target = "quality_percentile_modern",
                                       modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                       legacy_methods = c("ruleset1"),
                                       n_top_features = 20,
                                       xgb_params = NULL,
                                       quiet = FALSE) {

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required. Install with: install.packages('xgboost')")
  }

  if (!requireNamespace("shapviz", quietly = TRUE)) {
    stop("Package 'shapviz' is required. Install with: install.packages('shapviz')")
  }

  # ═══════════════════════════════════════════════════════════════
  # 1. Aggregate Data Across Genes
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Aggregating data across genes...")

  # Handle single analysis object
  if (inherits(batch_result, "mutateR_consensus_analysis")) {
    gene_name <- batch_result$metadata$gene_symbol %||% batch_result$metadata$gene_id %||% "gene"
    batch_result <- structure(
      list(batch_result),
      names = gene_name,
      class = c("mutateR_consensus_batch", "list")
    )
  }

  all_features <- list()
  all_targets <- list()

  for (gene in names(batch_result)) {
    obj <- batch_result[[gene]]
    if (is.null(obj$scores) || nrow(obj$scores) == 0) next

    # Extract features
    feat <- extract_grna_features(obj$scores)
    feat$gene <- gene

    # Compute targets with tiered structure
    targ <- compute_shap_targets(obj$scores,
                                 modern_methods = modern_methods,
                                 legacy_methods = legacy_methods)

    all_features[[gene]] <- feat
    all_targets[[gene]] <- targ
  }

  if (length(all_features) == 0) {
    stop("No valid score data found in batch_result")
  }

  features_df <- do.call(rbind, all_features)
  targets_df <- do.call(rbind, all_targets)

  if (!quiet) {
    message("Aggregated ", nrow(features_df), " gRNAs from ",
            length(all_features), " genes")
  }

  # ═══════════════════════════════════════════════════════════════
  # 2. Prepare Model Data
  # ═══════════════════════════════════════════════════════════════

  feature_cols <- setdiff(names(features_df), c("grna_id", "gene"))
  X <- as.matrix(features_df[, feature_cols])

  # Validate target exists
  if (!target %in% names(targets_df)) {
    available_targets <- names(targets_df)
    # Group targets for cleaner error message
    quality_targets <- grep("quality|score_", available_targets, value = TRUE)
    rank_targets <- grep("rank", available_targets, value = TRUE)
    other_targets <- setdiff(available_targets, c(quality_targets, rank_targets, "grna_id"))

    stop("Target '", target, "' not found.\n",
         "Quality targets: ", paste(head(quality_targets, 10), collapse = ", "), "\n",
         "Rank targets: ", paste(head(rank_targets, 10), collapse = ", "), "\n",
         "Other targets: ", paste(head(other_targets, 10), collapse = ", "))
  }

  y <- targets_df[[target]]

  # Remove rows with NA in features or target
  valid <- complete.cases(X) & !is.na(y)
  X <- X[valid, , drop = FALSE]
  y <- y[valid]
  features_df_valid <- features_df[valid, ]
  targets_df_valid <- targets_df[valid, ]

  if (sum(valid) < 100) {
    warning("Only ", sum(valid), " complete observations available. ",
            "Results may be unreliable.")
  }

  if (!quiet) message("Training on ", nrow(X), " complete observations")

  # ═══════════════════════════════════════════════════════════════
  # 3. Train XGBoost Model
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Training gradient boosting model...")

  dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

  # Default parameters (can be overridden)
  default_params <- list(
    objective = "reg:squarederror",
    max_depth = 6,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 5
  )

  if (!is.null(xgb_params)) {
    default_params <- modifyList(default_params, xgb_params)
  }

  # Train with early stopping via cross-validation to determine optimal rounds
  cv_result <- tryCatch({
    xgboost::xgb.cv(
      params = default_params,
      data = dtrain,
      nrounds = 500,
      nfold = 5,
      early_stopping_rounds = 20,
      verbose = 0
    )
  }, error = function(e) NULL)

  # Determine optimal number of rounds
  if (!is.null(cv_result)) {
    best_nrounds <- cv_result$best_iteration
    if (is.null(best_nrounds) || best_nrounds < 10) best_nrounds <- 100
  } else {
    best_nrounds <- 100
  }

  model <- xgboost::xgb.train(
    params = default_params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )

  # ═══════════════════════════════════════════════════════════════
  # 4. Compute SHAP Values
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Computing SHAP values...")

  shap_values <- shapviz::shapviz(model, X_pred = X, X = X)

  # ═══════════════════════════════════════════════════════════════
  # 5. Generate Visualisations
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Generating visualizations...")

  # Beeswarm plot (shows distribution of SHAP values)
  beeswarm_plot <- tryCatch({
    shapviz::sv_importance(shap_values, kind = "beeswarm",
                           max_display = n_top_features) +
      ggplot2::ggtitle(paste0("SHAP Feature Importance: ", target)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }, error = function(e) NULL)

  # Bar plot (shows mean |SHAP|)
  bar_plot <- tryCatch({
    shapviz::sv_importance(shap_values, kind = "bar",
                           max_display = n_top_features) +
      ggplot2::ggtitle(paste0("Mean |SHAP| Importance: ", target)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }, error = function(e) NULL)

  # ═══════════════════════════════════════════════════════════════
  # 6. Extract Feature Importance Table
  # ═══════════════════════════════════════════════════════════════

  # Get SHAP-based importance
  shap_matrix <- shapviz::get_shap_values(shap_values)
  mean_abs_shap <- colMeans(abs(shap_matrix))

  importance_df <- data.frame(
    feature = names(mean_abs_shap),
    mean_abs_shap = as.numeric(mean_abs_shap),
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(-importance_df$mean_abs_shap), ]
  importance_df$rank <- seq_len(nrow(importance_df))
  rownames(importance_df) <- NULL

  # Add direction information (mean SHAP value, not absolute)
  mean_shap <- colMeans(shap_matrix)
  importance_df$mean_shap <- mean_shap[importance_df$feature]
  importance_df$direction <- ifelse(importance_df$mean_shap > 0, "positive", "negative")

  # ═══════════════════════════════════════════════════════════════
  # 7. Compute Model Performance Metrics
  # ═══════════════════════════════════════════════════════════════

  predictions <- predict(model, X)

  # R-squared
  ss_res <- sum((y - predictions)^2)
  ss_tot <- sum((y - mean(y))^2)
  r_squared <- 1 - (ss_res / ss_tot)

  # RMSE
  rmse <- sqrt(mean((y - predictions)^2))

  # Correlation
  pred_cor <- cor(y, predictions, method = "spearman")

  # ═══════════════════════════════════════════════════════════════
  # 8. Return Results
  # ═══════════════════════════════════════════════════════════════

  result <- list(
    model = model,
    shap_values = shap_values,
    feature_importance = importance_df,
    plots = list(
      beeswarm = beeswarm_plot,
      bar = bar_plot
    ),
    performance = list(
      r_squared = r_squared,
      rmse = rmse,
      spearman_cor = pred_cor,
      n_rounds = best_nrounds
    ),
    data = list(
      features = features_df_valid,
      targets = targets_df_valid,
      X = X,
      y = y
    ),
    metadata = list(
      target = target,
      n_gRNAs = nrow(X),
      n_genes = length(all_features),
      n_features = ncol(X),
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      xgb_params = default_params
    )
  )

  class(result) <- c("mutateR_shap_analysis", class(result))

  if (!quiet) {
    message("Analysis complete. Model R² = ", round(r_squared, 3),
            ", Spearman ρ = ", round(pred_cor, 3))
  }

  return(result)
}


#' Print method for SHAP analysis results
#' @export
print.mutateR_shap_analysis <- function(x, ...) {

  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║         mutateR SHAP Feature Importance Analysis              ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

  cat("TARGET VARIABLE\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("Target:          ", x$metadata$target, "\n")

  # Add interpretation based on target type
  target_interp <- get_target_interpretation(x$metadata$target)
  if (!is.null(target_interp)) {
    cat("Interpretation:  ", target_interp, "\n")
  }
  cat("\n")

  cat("DATASET\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("gRNAs analysed:  ", x$metadata$n_gRNAs, "\n")
  cat("Genes included:  ", x$metadata$n_genes, "\n")
  cat("Features:        ", x$metadata$n_features, "\n")
  cat("Modern methods:  ", paste(x$metadata$modern_methods, collapse = ", "), "\n")
  cat("Legacy methods:  ", paste(x$metadata$legacy_methods, collapse = ", "), "\n")
  cat("\n")

  cat("MODEL PERFORMANCE\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("R²:              %.3f\n", x$performance$r_squared))
  cat(sprintf("RMSE:            %.3f\n", x$performance$rmse))
  cat(sprintf("Spearman ρ:      %.3f\n", x$performance$spearman_cor))
  cat(sprintf("XGBoost rounds:  %d\n", x$performance$n_rounds))
  cat("\n")

  cat("TOP 15 FEATURES (by mean |SHAP|)\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  top15 <- head(x$feature_importance, 15)
  top15$mean_abs_shap <- sprintf("%.4f", top15$mean_abs_shap)
  top15$mean_shap <- sprintf("%+.4f", top15$mean_shap)
  print(top15[, c("rank", "feature", "mean_abs_shap", "direction")], row.names = FALSE)

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Access plots:           $plots$beeswarm, $plots$bar\n")
  cat("Access SHAP values:     $shap_values\n")
  cat("Access full importance: $feature_importance\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  invisible(x)
}


#' Get human-readable interpretation for target variables
#' @noRd
get_target_interpretation <- function(target) {
  interpretations <- list(
    quality_percentile_modern = "Higher = better efficacy (modern method consensus)",
    quality_percentile_robust = "Higher = better efficacy (outlier-resistant)",
    quality_percentile_all = "Higher = better efficacy (all methods)",
    quality_inverted_rank_modern = "Higher = better efficacy (modern methods)",
    concordance_modern = "Higher = more agreement among modern methods",
    discordance_modern = "Higher = more disagreement among modern methods",
    rank_variance_modern = "Higher = more disagreement among modern methods",
    ruleset1_divergence = "Positive = RuleSet1 ranks worse than modern consensus",
    ruleset1_abs_divergence = "Higher = larger RuleSet1 vs modern disagreement",
    mean_rank_modern = "Lower = better efficacy (modern methods)",
    mean_rank_all = "Lower = better efficacy (all methods)"
  )

  # Check for exact match
  if (target %in% names(interpretations)) {
    return(interpretations[[target]])
  }

  # Check for method-specific scores
  if (grepl("^score_", target)) {
    method <- sub("^score_", "", target)
    return(paste0("Higher = better efficacy per ", method))
  }

  if (grepl("^rank_", target)) {
    method <- sub("^rank_", "", target)
    return(paste0("Lower = better efficacy per ", method))
  }

  return(NULL)
}


#' Generate SHAP dependence plots for specific features
#'
#' Creates dependence plots showing the relationship between a feature's
#' value and its SHAP contribution, optionally colored by an interaction feature.
#'
#' @param shap_result Output from analyze_feature_importance()
#' @param features Character vector of feature names to plot
#' @param color_by Character. Feature to use for coloring points (default: auto-select)
#' @param ncol Integer. Number of columns in plot grid (default: 2)
#'
#' @return A ggplot object or list of ggplot objects
#' @export
plot_shap_dependence <- function(shap_result,
                                 features = NULL,
                                 color_by = "auto",
                                 ncol = 2) {

  if (!inherits(shap_result, "mutateR_shap_analysis")) {
    stop("Input must be a mutateR_shap_analysis object")
  }

  # Default to top 4 features
  if (is.null(features)) {
    features <- head(shap_result$feature_importance$feature, 4)
  }

  # Validate features exist
  available <- colnames(shap_result$data$X)
  invalid <- setdiff(features, available)
  if (length(invalid) > 0) {
    stop("Features not found: ", paste(invalid, collapse = ", "))
  }

  # Generate plots
  plots <- lapply(features, function(feat) {
    tryCatch({
      shapviz::sv_dependence(shap_result$shap_values, v = feat, color_var = color_by) +
        ggplot2::ggtitle(feat)
    }, error = function(e) {
      message("Failed to create dependence plot for ", feat, ": ", e$message)
      NULL
    })
  })

  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) {
    warning("No dependence plots could be generated")
    return(NULL)
  }

  if (length(plots) == 1) {
    return(plots[[1]])
  }

  # Combine with patchwork if available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(plots, ncol = ncol))
  }

  return(plots)
}


#' Compare SHAP results between two analyses
#'
#' Generates comparison visualisations and statistics for two SHAP analyses,
#' useful for comparing feature importance between targets (e.g., modern quality
#' vs. legacy divergence) or between different datasets.
#'
#' @param shap1 First mutateR_shap_analysis object
#' @param shap2 Second mutateR_shap_analysis object
#' @param label1 Character. Label for first analysis (default: uses target name)
#' @param label2 Character. Label for second analysis (default: uses target name)
#' @param n_features Integer. Number of top features to compare (default: 20)
#'
#' @return List containing comparison data frame and visualisation
#' @export
compare_shap_analyses <- function(shap1, shap2,
                                  label1 = NULL, label2 = NULL,
                                  n_features = 20) {

  if (is.null(label1)) label1 <- shap1$metadata$target
  if (is.null(label2)) label2 <- shap2$metadata$target

  # Extract importance rankings
  imp1 <- shap1$feature_importance[, c("feature", "mean_abs_shap", "rank")]
  imp2 <- shap2$feature_importance[, c("feature", "mean_abs_shap", "rank")]

  names(imp1)[2:3] <- paste0(c("importance_", "rank_"), label1)
  names(imp2)[2:3] <- paste0(c("importance_", "rank_"), label2)

  # Merge
  comparison <- merge(imp1, imp2, by = "feature", all = TRUE)

  # Calculate rank difference
  rank_col1 <- paste0("rank_", label1)
  rank_col2 <- paste0("rank_", label2)

  comparison$rank_diff <- comparison[[rank_col1]] - comparison[[rank_col2]]
  comparison$abs_rank_diff <- abs(comparison$rank_diff)

  # Classify features
  comparison$category <- sapply(1:nrow(comparison), function(i) {
    r1 <- comparison[[rank_col1]][i]
    r2 <- comparison[[rank_col2]][i]

    top_threshold <- n_features

    if (is.na(r1) && !is.na(r2) && r2 <= top_threshold) {
      return(paste0("Unique to ", label2))
    }
    if (!is.na(r1) && is.na(r2) && r1 <= top_threshold) {
      return(paste0("Unique to ", label1))
    }
    if (!is.na(r1) && !is.na(r2)) {
      if (r1 <= top_threshold && r2 <= top_threshold) {
        return("Important in both")
      }
      if (r1 <= top_threshold) {
        return(paste0("More important for ", label1))
      }
      if (r2 <= top_threshold) {
        return(paste0("More important for ", label2))
      }
    }
    return("Less important in both")
  })

  # Sort by minimum rank (most important features first)
  comparison$min_rank <- pmin(comparison[[rank_col1]], comparison[[rank_col2]], na.rm = TRUE)
  comparison <- comparison[order(comparison$min_rank), ]
  rownames(comparison) <- NULL

  # Create visualisation
  plot_df <- head(comparison, n_features * 2)
  plot_df <- plot_df[plot_df$category != "Less important in both", ]

  if (nrow(plot_df) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {

    # Rank comparison scatter plot
    rank_plot <- ggplot2::ggplot(plot_df,
                                 ggplot2::aes_string(x = rank_col1, y = rank_col2,
                                                     color = "category")) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_text(ggplot2::aes(label = feature),
                         size = 2.5, hjust = -0.1, vjust = -0.1,
                         check_overlap = TRUE) +
      ggplot2::scale_color_brewer(palette = "Set2") +
      ggplot2::labs(title = paste("Feature Rank Comparison:", label1, "vs", label2),
                    x = paste("Rank in", label1),
                    y = paste("Rank in", label2),
                    color = "Category") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

  } else {
    rank_plot <- NULL
  }

  return(list(
    comparison = comparison,
    plot = rank_plot,
    summary = list(
      n_shared_top = sum(comparison$category == "Important in both", na.rm = TRUE),
      correlation = cor(comparison[[rank_col1]], comparison[[rank_col2]],
                        use = "complete.obs", method = "spearman")
    )
  ))
}
