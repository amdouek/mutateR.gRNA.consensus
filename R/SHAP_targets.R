#' Compute target variables for tiered SHAP analysis
#'
#' Calculates outcome variables with explicit handling of modern vs. legacy
#' methods. Provides tiered consensus metrics for flexible analysis that
#' prevents legacy methods from dominating variance-based statistics.
#'
#' @param scores_df Data frame from score_grnas_multi()
#' @param modern_methods Character vector. Modern scoring methods.
#'        Default: c("deepspcas9", "deephf", "ruleset3") for Cas9,
#'                 c("deepcpf1", "enpamgb") for Cas12a.
#' @param legacy_methods Character vector. Legacy methods to model separately.
#'        Default: c("ruleset1") for Cas9, empty for Cas12a.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{grna_id}{gRNA identifier}
#'     \item{mean_score_modern}{Mean normalized score (modern methods)}
#'     \item{mean_rank_modern}{Mean rank across modern methods (1 = best)}
#'     \item{max_rank_modern}{Worst rank among modern methods}
#'     \item{min_rank_modern}{Best rank among modern methods}
#'     \item{rank_variance_modern}{Variance in ranks (modern methods)}
#'     \item{rank_range_modern}{Range of ranks (modern methods)}
#'     \item{quality_percentile_modern}{Percentile score (higher = better)}
#'     \item{concordance_modern}{Agreement score (higher = more agreement)}
#'     \item{discordance_modern}{Disagreement score (higher = more disagreement)}
#'     \item{median_rank_all}{Median rank across all methods}
#'     \item{trimmed_mean_rank_all}{Trimmed mean rank (robust)}
#'     \item{quality_percentile_robust}{Median-based quality percentile}
#'     \item{<legacy>_divergence}{Legacy method divergence from modern consensus}
#'     \item{is_top20/50/100_modern}{Binary flags for top-ranked gRNAs}
#'     \item{score_<method>}{Raw scores for each method}
#'     \item{rank_<method>}{Ranks for each method}
#'   }
#'
#' @details
#' The tiered approach addresses a key challenge: legacy methods (like RuleSet1)
#' may systematically differ from modern deep learning methods, which can
#' dominate simple variance-based consensus metrics.
#'
#' \strong{Tier 1 - Modern Method Consensus:}
#' Metrics computed only among modern methods (DeepSpCas9, DeepHF, RuleSet3).
#' These methods show higher inter-correlation (~0.70) and likely better
#' reflect true gRNA efficacy.
#'
#' \strong{Tier 2 - Robust Full Consensus:}
#' Metrics using all methods but with robust statistics (median, trimmed mean)
#' that resist influence from a single outlier method.
#'
#' \strong{Tier 3 - Legacy Divergence:}
#' Explicit modeling of where legacy methods disagree with modern consensus.
#' Positive divergence means legacy ranks the gRNA worse than modern methods.
#'
#' \strong{Tier 4 - Binary Selection Flags:}
#' Useful for classification tasks (e.g., "is this a top-50 gRNA?").
#'
#' @examples
#' \dontrun{
#' # Get multi-method scores
#' scores <- score_grnas_multi(sites)
#'
#' # Compute all target variables
#' targets <- compute_shap_targets(scores)
#'
#' # View available targets
#' names(targets)
#'
#' # Use modern quality for primary analysis
#' y <- targets$quality_percentile_modern
#'
#' # Analyse legacy divergence separately
#' y_legacy <- targets$ruleset1_divergence
#' }
#'
#' @seealso
#' \code{\link{extract_grna_features}} for feature computation,
#' \code{\link{analyze_feature_importance}} for SHAP analysis,
#' \code{\link{get_target_interpretation}} for target descriptions
#'
#' @export
compute_shap_targets <- function(scores_df,
                                 modern_methods = NULL,
                                 legacy_methods = NULL) {

  # Identify Available Methods
  all_available <- attr(scores_df, "methods")
  nuclease <- attr(scores_df, "nuclease")

  # Fallback: detect numeric columns that aren't metadata
  if (is.null(all_available)) {
    meta_cols <- c("grna_id", "protospacer_sequence", "sequence_context",
                   "exon_rank", "chr", "start", "end", "strand")
    all_available <- setdiff(names(scores_df), meta_cols)
    all_available <- all_available[vapply(scores_df[all_available], is.numeric, logical(1))]
  }

  # Set Default Method Classification
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

  # Initialize output
  targets <- data.frame(
    grna_id = if ("grna_id" %in% names(scores_df)) {
      scores_df$grna_id
    } else {
      paste0("gRNA_", seq_len(n_grnas))
    },
    stringsAsFactors = FALSE
  )

  # TIER 1: Modern Method Consensus
  if (length(modern_methods) >= 2) {
    modern_scores <- scores_df[, modern_methods, drop = FALSE]

    # Normalize scores to 0-1 range
    modern_norm <- as.data.frame(lapply(modern_scores, function(x) {
      rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      if (rng == 0) return(rep(0.5, length(x)))
      (x - min(x, na.rm = TRUE)) / rng
    }))

    # Compute ranks (1 = best/highest score)
    modern_ranks <- as.data.frame(lapply(modern_scores, function(x) {
      rank(-x, na.last = "keep", ties.method = "average")
    }))

    # Summary statistics
    targets$mean_score_modern <- rowMeans(modern_norm, na.rm = TRUE)
    targets$mean_rank_modern <- rowMeans(modern_ranks, na.rm = TRUE)
    targets$max_rank_modern <- apply(modern_ranks, 1, max, na.rm = TRUE)
    targets$min_rank_modern <- apply(modern_ranks, 1, min, na.rm = TRUE)
    targets$rank_variance_modern <- apply(modern_ranks, 1, var, na.rm = TRUE)
    targets$rank_range_modern <- targets$max_rank_modern - targets$min_rank_modern

    # Quality percentile (higher = better, 0-100 scale)
    targets$quality_percentile_modern <- 100 * (1 - (targets$mean_rank_modern - 1) / (n_grnas - 1))

    # Inverted rank for intuitive interpretation (higher = better)
    targets$quality_inverted_rank_modern <- (n_grnas + 1) - targets$mean_rank_modern

    # Concordance score (higher = more agreement, 0-100 scale)
    max_var <- max(targets$rank_variance_modern, na.rm = TRUE)
    targets$concordance_modern <- if (max_var > 0) {
      100 * (1 - targets$rank_variance_modern / max_var)
    } else {
      rep(100, n_grnas)
    }

    # Discordance score (higher = more disagreement, 0-100 scale)
    targets$discordance_modern <- if (max_var > 0) {
      100 * (targets$rank_variance_modern / max_var)
    } else {
      rep(0, n_grnas)
    }
  }

  # TIER 2: Robust Full Consensus (All Methods)
  all_scores <- scores_df[, all_methods, drop = FALSE]

  # Normalize all scores
  all_norm <- as.data.frame(lapply(all_scores, function(x) {
    rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / rng
  }))

  # Compute ranks for all methods
  all_ranks <- as.data.frame(lapply(all_scores, function(x) {
    rank(-x, na.last = "keep", ties.method = "average")
  }))

  # Robust statistics (resistant to single outlier method)
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

  # TIER 3: Legacy Divergence
  if (length(legacy_methods) >= 1 && length(modern_methods) >= 2) {
    for (legacy in legacy_methods) {
      legacy_rank <- all_ranks[[legacy]]
      col_prefix <- gsub("[^a-z0-9]", "_", tolower(legacy))

      # Signed divergence: positive = legacy ranks it worse than modern consensus
      targets[[paste0(col_prefix, "_divergence")]] <-
        legacy_rank - targets$mean_rank_modern

      # Absolute divergence (magnitude regardless of direction)
      targets[[paste0(col_prefix, "_abs_divergence")]] <-
        abs(targets[[paste0(col_prefix, "_divergence")]])

      # Scaled divergence (percentage of total gRNAs)
      targets[[paste0(col_prefix, "_divergence_pct")]] <-
        100 * targets[[paste0(col_prefix, "_divergence")]] / n_grnas

      # Direction flag: is legacy more optimistic (ranks it better)?
      targets[[paste0(col_prefix, "_more_optimistic")]] <-
        as.integer(legacy_rank < targets$mean_rank_modern)
    }
  }

  # TIER 4: Binary Selection Flags
  # Modern method top-N flags
  if ("max_rank_modern" %in% names(targets)) {
    targets$is_top20_modern <- as.integer(targets$max_rank_modern <= 20)
    targets$is_top50_modern <- as.integer(targets$max_rank_modern <= 50)
    targets$is_top100_modern <- as.integer(targets$max_rank_modern <= 100)
  }

  # Robust top-N flag
  targets$is_top50_robust <- as.integer(targets$median_rank_all <= 50)

  # High discordance flag (top 10% most discordant)
  if ("rank_variance_modern" %in% names(targets)) {
    targets$is_discordant_modern <- as.integer(
      targets$rank_variance_modern > quantile(targets$rank_variance_modern, 0.9, na.rm = TRUE)
    )
  }

  # Legacy divergent flag (top 10% most divergent from modern)
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

  # Individual Method Scores and Ranks
  for (m in all_methods) {
    targets[[paste0("score_", m)]] <- all_scores[[m]]
    targets[[paste0("rank_", m)]] <- all_ranks[[m]]
  }

  # Add Metadata as Attributes
  attr(targets, "modern_methods") <- modern_methods
  attr(targets, "legacy_methods") <- legacy_methods
  attr(targets, "n_grnas") <- n_grnas
  attr(targets, "nuclease") <- nuclease

  return(targets)
}


#' Get human-readable interpretation for target variables
#'
#' Returns a plain-language description of what a target variable represents
#' and how to interpret its values (higher/lower = better/worse).
#'
#' @param target Character. Name of target variable.
#'
#' @return Character string with interpretation, or NULL if target not recognized.
#'
#' @examples
#' get_target_interpretation("quality_percentile_modern")
#' # "Higher = better efficacy (modern method consensus)"
#'
#' get_target_interpretation("ruleset1_divergence")
#' # "Positive = RuleSet1 ranks worse than modern consensus"
#'
#' get_target_interpretation("score_deepspcas9")
#' # "Higher = better efficacy per deepspcas9"
#'
#' @export
get_target_interpretation <- function(target) {

  # Static interpretations for known targets
  interpretations <- list(
    # Quality metrics (higher = better)
    quality_percentile_modern = "Higher = better efficacy (modern method consensus)",
    quality_percentile_robust = "Higher = better efficacy (outlier-resistant)",
    quality_percentile_all = "Higher = better efficacy (all methods)",
    quality_inverted_rank_modern = "Higher = better efficacy (modern methods)",
    mean_score_modern = "Higher = better efficacy (modern methods)",

    # Rank metrics (lower = better)
    mean_rank_modern = "Lower = better efficacy (modern methods)",
    mean_rank_all = "Lower = better efficacy (all methods)",
    median_rank_all = "Lower = better efficacy (median across all)",
    trimmed_mean_rank_all = "Lower = better efficacy (trimmed mean)",
    max_rank_modern = "Lower = more consistent top performer (modern)",
    max_rank_all = "Lower = more consistent top performer (all)",

    # Agreement metrics
    concordance_modern = "Higher = more agreement among modern methods",
    discordance_modern = "Higher = more disagreement among modern methods",
    rank_variance_modern = "Higher = more disagreement among modern methods",
    rank_variance_all = "Higher = more disagreement among all methods",
    rank_range_modern = "Higher = larger spread in rankings (modern)",
    rank_iqr_all = "Higher = larger spread in rankings (all)",

    # Legacy divergence
    ruleset1_divergence = "Positive = RuleSet1 ranks worse than modern consensus",
    ruleset1_abs_divergence = "Higher = larger RuleSet1 vs modern disagreement",
    ruleset1_divergence_pct = "Positive % = RuleSet1 relatively pessimistic",
    ruleset1_more_optimistic = "1 = RuleSet1 ranks better than modern consensus",

    # Binary flags
    is_top20_modern = "1 = ranked in top 20 by all modern methods",
    is_top50_modern = "1 = ranked in top 50 by all modern methods",
    is_top100_modern = "1 = ranked in top 100 by all modern methods",
    is_top50_robust = "1 = median rank across all methods <= 50",
    is_discordant_modern = "1 = top 10% most disagreed upon (modern)",
    is_legacy_divergent = "1 = top 10% largest legacy divergence"
  )

  # Check for exact match
  if (target %in% names(interpretations)) {
    return(interpretations[[target]])
  }

  # Check for method-specific scores (score_<method>)
  if (grepl("^score_", target)) {
    method <- sub("^score_", "", target)
    return(paste0("Higher = better efficacy per ", method))
  }

  # Check for method-specific ranks (rank_<method>)
  if (grepl("^rank_", target)) {
    method <- sub("^rank_", "", target)
    return(paste0("Lower = better efficacy per ", method))
  }

  # Unknown target
  return(NULL)
}


#' List all available SHAP target variables
#'
#' Returns a data frame describing all target variables that can be
#' computed by \code{compute_shap_targets()}.
#'
#' @param include_method_specific Logical. Include score_* and rank_* columns
#'        for individual methods (default FALSE).
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{target}{Target variable name}
#'     \item{tier}{Tier classification (1-4)}
#'     \item{category}{Functional category}
#'     \item{interpretation}{Plain-language description}
#'   }
#'
#' @examples
#' targets <- list_shap_targets()
#' print(targets)
#'
#' # Filter to quality metrics
#' quality_targets <- targets[targets$category == "quality", ]
#'
#' @export
list_shap_targets <- function(include_method_specific = FALSE) {

  targets <- data.frame(
    target = character(),
    tier = integer(),
    category = character(),
    interpretation = character(),
    stringsAsFactors = FALSE
  )

  # Tier 1: Modern Method Consensus
  tier1 <- data.frame(
    target = c("mean_score_modern", "mean_rank_modern", "max_rank_modern",
               "min_rank_modern", "rank_variance_modern", "rank_range_modern",
               "quality_percentile_modern", "quality_inverted_rank_modern",
               "concordance_modern", "discordance_modern"),
    tier = 1L,
    category = c("quality", "quality", "quality", "quality",
                 "agreement", "agreement", "quality", "quality",
                 "agreement", "agreement"),
    stringsAsFactors = FALSE
  )

  # Tier 2: Robust Full Consensus
  tier2 <- data.frame(
    target = c("median_rank_all", "trimmed_mean_rank_all", "rank_iqr_all",
               "max_rank_all", "mean_rank_all", "rank_variance_all",
               "quality_percentile_all", "quality_percentile_robust"),
    tier = 2L,
    category = c("quality", "quality", "agreement", "quality",
                 "quality", "agreement", "quality", "quality"),
    stringsAsFactors = FALSE
  )

  # Tier 3: Legacy Divergence
  tier3 <- data.frame(
    target = c("ruleset1_divergence", "ruleset1_abs_divergence",
               "ruleset1_divergence_pct", "ruleset1_more_optimistic"),
    tier = 3L,
    category = c("divergence", "divergence", "divergence", "divergence"),
    stringsAsFactors = FALSE
  )

  # Tier 4: Binary Flags
  tier4 <- data.frame(
    target = c("is_top20_modern", "is_top50_modern", "is_top100_modern",
               "is_top50_robust", "is_discordant_modern", "is_legacy_divergent"),
    tier = 4L,
    category = c("selection", "selection", "selection",
                 "selection", "selection", "selection"),
    stringsAsFactors = FALSE
  )

  targets <- rbind(tier1, tier2, tier3, tier4)

  # Add interpretations
  targets$interpretation <- vapply(targets$target, function(t) {
    interp <- get_target_interpretation(t)
    if (is.null(interp)) "" else interp
  }, character(1))

  # Add method-specific if requested
  if (include_method_specific) {
    methods <- c("ruleset1", "deepspcas9", "ruleset3", "deephf")

    method_targets <- data.frame(
      target = c(paste0("score_", methods), paste0("rank_", methods)),
      tier = NA_integer_,
      category = rep(c("score", "rank"), each = length(methods)),
      stringsAsFactors = FALSE
    )

    method_targets$interpretation <- vapply(method_targets$target, function(t) {
      interp <- get_target_interpretation(t)
      if (is.null(interp)) "" else interp
    }, character(1))

    targets <- rbind(targets, method_targets)
  }

  return(targets)
}


#' Recommend target variable for analysis
#'
#' Suggests the most appropriate target variable based on the analysis goal.
#'
#' @param goal Character. Analysis goal, one of:
#'   \describe{
#'     \item{"quality"}{Primary efficacy prediction (default)}
#'     \item{"agreement"}{Method disagreement analysis}
#'     \item{"legacy"}{Legacy method divergence}
#'     \item{"robust"}{Outlier-resistant quality}
#'     \item{"selection"}{Binary top-gRNA classification}
#'   }
#'
#' @return Character. Recommended target variable name.
#'
#' @examples
#' recommend_shap_target("quality")
#' # "quality_percentile_modern"
#'
#' recommend_shap_target("legacy")
#' # "ruleset1_divergence"
#'
#' @export
recommend_shap_target <- function(goal = c("quality", "agreement", "legacy",
                                           "robust", "selection")) {

  goal <- match.arg(goal)

  recommendations <- list(
    quality = "quality_percentile_modern",
    agreement = "discordance_modern",
    legacy = "ruleset1_divergence",
    robust = "quality_percentile_robust",
    selection = "is_top50_modern"
  )

  target <- recommendations[[goal]]

  message("Recommended target: ", target)
  message("Interpretation: ", get_target_interpretation(target))

  return(target)
}
