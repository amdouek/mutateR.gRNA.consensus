#' @title SHAP Analysis for gRNA Scoring Methods
#' @description
#' Core SHAP analysis functions for identifying sequence features that drive
#' gRNA scoring method predictions and disagreements. Uses XGBoost gradient
#' boosting with SHAP (SHapley Additive exPlanations) values.
#'
#' @name SHAP_analysis
NULL


#' Run SHAP analysis for gRNA features
#'
#' Uses a gradient boosting model with SHAP values to identify which sequence
#' features drive high scores and cross-method agreement/disagreement.
#' Supports tiered target variables for nuanced analysis.
#'
#' @param batch_result A mutateR_consensus_batch or mutateR_consensus_analysis object
#' @param target Character. Target variable to model. Common options:
#'   \describe{
#'     \item{"quality_percentile_modern"}{Modern method consensus (recommended primary)}
#'     \item{"quality_percentile_robust"}{Outlier-resistant quality}
#'     \item{"discordance_modern"}{Disagreement across modern methods}
#'     \item{"ruleset1_divergence"}{Legacy method divergence}
#'     \item{"score_<method>"}{Individual method score}
#'   }
#' @param modern_methods Character vector of modern scoring methods.
#'        Default: c("deepspcas9", "deephf", "ruleset3")
#' @param legacy_methods Character vector of legacy scoring methods.
#'        Default: c("ruleset1")
#' @param n_top_features Integer. Number of top features to display (default 20)
#' @param xgb_params List. XGBoost parameters (optional, uses sensible defaults)
#' @param rasterize Logical. Rasterise beeswarm points for large datasets.
#'        Default NULL means auto-determine (TRUE if n > 5000).
#' @param raster_dpi Integer. DPI for rasterized layers (default 300).
#' @param subsample_beeswarm Integer or NULL. If set, subsample points for beeswarm plot.
#'        NULL (default) means auto-determine based on dataset size.
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return List of class "mutateR_shap_analysis" containing:
#'   \describe{
#'     \item{model}{Trained XGBoost model}
#'     \item{shap_values}{shapviz object with SHAP values}
#'     \item{feature_importance}{Data frame of features ranked by importance}
#'     \item{plots}{List of ggplot objects (beeswarm, bar, waterfall)}
#'     \item{performance}{List with R², RMSE, Spearman correlation}
#'     \item{data}{List with features, targets, X matrix, y vector}
#'     \item{metadata}{List with analysis parameters}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Aggregates gRNA data across all genes in batch_result
#'   \item Extracts sequence features using \code{\link{extract_grna_features}}
#'   \item Computes target variables using \code{\link{compute_shap_targets}}
#'   \item Trains XGBoost model with cross-validation for early stopping
#'   \item Computes SHAP values for feature importance
#'   \item Generates visualization plots
#' }
#'
#' XGBoost default parameters are tuned for gRNA data:
#' \itemize{
#'   \item max_depth = 6 (moderate tree depth)
#'   \item eta = 0.1 (learning rate)
#'   \item subsample = 0.8 (row sampling)
#'   \item colsample_bytree = 0.8 (feature sampling)
#'   \item min_child_weight = 5 (minimum leaf size)
#' }
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
#' # SHAP analysis for quality
#' shap_quality <- analyze_feature_importance(
#'   batch,
#'   target = "quality_percentile_modern"
#' )
#'
#' # View results
#' print(shap_quality)
#' shap_quality$plots$beeswarm
#'
#' # SHAP analysis for disagreement
#' shap_discord <- analyze_feature_importance(
#'   batch,
#'   target = "discordance_modern"
#' )
#' }
#'
#' @seealso
#' \code{\link{run_tiered_shap_analysis}} for multi-tier analysis,
#' \code{\link{extract_grna_features}} for feature computation,
#' \code{\link{compute_shap_targets}} for target variables
#'
#' @export
analyze_feature_importance <- function(batch_result,
                                       target = "quality_percentile_modern",
                                       modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                       legacy_methods = c("ruleset1"),
                                       n_top_features = 20,
                                       xgb_params = NULL,
                                       rasterize = NULL,
                                       raster_dpi = 300,
                                       subsample_beeswarm = NULL,
                                       quiet = FALSE) {

  # Check Dependencies
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required. Install with: install.packages('xgboost')")
  }

  if (!requireNamespace("shapviz", quietly = TRUE)) {
    stop("Package 'shapviz' is required. Install with: install.packages('shapviz')")
  }

  # 1. Aggregate Data Across Genes
  if (!quiet) message("Aggregating data across genes...")

  # Handle single analysis object
  if (inherits(batch_result, "mutateR_consensus_analysis")) {
    gene_name <- batch_result$metadata$gene_symbol %||%
      batch_result$metadata$gene_id %||%
      "gene"
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
    targ <- compute_shap_targets(
      obj$scores,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods
    )

    all_features[[gene]] <- feat
    all_targets[[gene]] <- targ
  }

  if (length(all_features) == 0) {
    stop("No valid score data found in batch_result")
  }

  features_df <- do.call(rbind, all_features)
  targets_df <- do.call(rbind, all_targets)
  n_genes <- length(all_features)

  if (!quiet) {
    message("Aggregated ", format(nrow(features_df), big.mark = ","),
            " gRNAs from ", n_genes, " genes")
  }

  # 2. Prepare Model Data
  feature_cols <- setdiff(names(features_df), c("grna_id", "gene"))
  X <- as.matrix(features_df[, feature_cols])

  # Validate target exists
  if (!target %in% names(targets_df)) {
    available_targets <- names(targets_df)
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

  n_gRNAs <- nrow(X)

  if (n_gRNAs < 100) {
    warning("Only ", n_gRNAs, " complete observations available. ",
            "Results may be unreliable.")
  }

  if (!quiet) {
    message("Training on ", format(n_gRNAs, big.mark = ","), " complete observations")
  }

  # 3. Train XGBoost Model
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

  # Train with early stopping via cross-validation
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

  # 4. Compute SHAP Values
  if (!quiet) message("Computing SHAP values...")

  shap_values <- shapviz::shapviz(model, X_pred = X, X = X)

  # 5. Compute Model Performance Metrics
  predictions <- predict(model, X)

  # R-squared
  ss_res <- sum((y - predictions)^2)
  ss_tot <- sum((y - mean(y))^2)
  r_squared <- 1 - (ss_res / ss_tot)

  # RMSE
  rmse <- sqrt(mean((y - predictions)^2))

  # Spearman correlation
  pred_cor <- cor(y, predictions, method = "spearman")

  # 6. Generate Visualisations
  if (!quiet) message("Generating visualizations...")

  plots <- generate_shap_plots(
    shap_values = shap_values,
    target = target,
    n_gRNAs = n_gRNAs,
    n_genes = n_genes,
    r_squared = r_squared,
    spearman_cor = pred_cor,
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    rasterize = rasterize,
    raster_dpi = raster_dpi,
    subsample_beeswarm = subsample_beeswarm
  )

  # 7. Extract Feature Importance Table
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

  # Add direction information
  mean_shap <- colMeans(shap_matrix)
  importance_df$mean_shap <- mean_shap[importance_df$feature]
  importance_df$direction <- ifelse(importance_df$mean_shap > 0, "positive", "negative")

  # 8. Return Results
  result <- list(
    model = model,
    shap_values = shap_values,
    feature_importance = importance_df,
    plots = plots,
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
      n_gRNAs = n_gRNAs,
      n_genes = n_genes,
      n_features = ncol(X),
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      xgb_params = default_params
    )
  )

  class(result) <- c("mutateR_shap_analysis", class(result))

  if (!quiet) {
    message("Analysis complete. Model R\u00B2 = ", round(r_squared, 3),
            ", Spearman \u03C1 = ", round(pred_cor, 3))
  }

  return(result)
}


#' Run tiered SHAP analysis across consensus metrics
#'
#' Executes SHAP analysis for multiple target variables representing different
#' aspects of gRNA quality and method agreement.
#'
#' @param batch_result A mutateR_consensus_batch object
#' @param modern_methods Character vector. Modern scoring methods.
#'        Default: c("deepspcas9", "deephf", "ruleset3")
#' @param legacy_methods Character vector. Legacy methods.
#'        Default: c("ruleset1")
#' @param n_top_features Integer. Number of top features to extract (default 25)
#' @param include_method_specific Logical. Also run SHAP for individual methods (default FALSE)
#' @param rasterize Logical. Rasterize beeswarm plots (default TRUE if n > 5000)
#' @param raster_dpi Integer. DPI for rasterisation (default 300)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return A list of class "mutateR_tiered_shap" containing:
#'   \describe{
#'     \item{modern_quality}{SHAP for quality_percentile_modern}
#'     \item{modern_discord}{SHAP for discordance_modern}
#'     \item{legacy_divergence}{SHAP for ruleset1_divergence}
#'     \item{robust_quality}{SHAP for quality_percentile_robust}
#'     \item{comparison}{Feature ranking comparison across analyses}
#'     \item{method_specific}{(optional) Individual method SHAP results}
#'     \item{metadata}{Analysis parameters and completion status}
#'   }
#'
#' @details
#' The four tiers address different analytical questions:
#'
#' \strong{modern_quality:} What features predict high efficacy according to
#' modern deep learning methods?
#'
#' \strong{modern_discord:} What features cause modern methods to disagree?
#'
#' \strong{legacy_divergence:} What features does legacy methods "miss" that modern
#' methods capture?
#'
#' \strong{robust_quality:} What features predict quality using outlier-resistant
#' statistics? Provides a conservative consensus view.
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Run batch analysis
#' batch <- analyze_score_consensus(
#'   c("TP53", "BRCA1", "EGFR", "MYC", "KRAS"),
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
#' # Compare feature rankings across tiers
#' shap_results$comparison
#' }
#'
#' @seealso
#' \code{\link{analyze_feature_importance}} for single-target analysis,
#' \code{\link{compare_shap_analyses}} for pairwise comparison
#'
#' @export
run_tiered_shap_analysis <- function(batch_result,
                                     modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                     legacy_methods = c("ruleset1"),
                                     n_top_features = 25,
                                     include_method_specific = FALSE,
                                     rasterize = NULL,
                                     raster_dpi = 300,
                                     quiet = FALSE) {

  # Validate input
  if (!inherits(batch_result, "mutateR_consensus_batch") &&
      !inherits(batch_result, "mutateR_consensus_analysis")) {
    stop("Input must be a mutateR_consensus_batch or mutateR_consensus_analysis object")
  }

  # Convert single analysis to batch format
  if (inherits(batch_result, "mutateR_consensus_analysis")) {
    gene_name <- batch_result$metadata$gene_symbol %||%
      batch_result$metadata$gene_id
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
    cat("\u2554", paste(rep("\u2550", 63), collapse = ""), "\u2557\n", sep = "")
    cat("\u2551           mutateR Tiered SHAP Analysis                        \u2551\n")
    cat("\u255A", paste(rep("\u2550", 63), collapse = ""), "\u255D\n\n", sep = "")
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
        n_top_features = n_top_features,
        rasterize = rasterize,
        raster_dpi = raster_dpi,
        quiet = TRUE
      )
      res$interpretation <- tier_info$interpretation
      res$tier_name <- tier_name

      if (!quiet) cat("    \u2713 Complete\n\n")
      res
    }, error = function(e) {
      if (!quiet) cat(sprintf("    \u2717 Failed: %s\n\n", e$message))
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
          n_top_features = n_top_features,
          rasterize = rasterize,
          raster_dpi = raster_dpi,
          quiet = TRUE
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
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")
    cat("Tiered SHAP analysis complete.\n")
    cat("Access results via: $modern_quality, $modern_discord,\n")
    cat("                    $legacy_divergence, $robust_quality\n")
    cat("Compare features via: $comparison\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")
  }

  return(results)
}


#' Compare feature rankings across SHAP tiers
#'
#' Internal function that compares feature importance rankings across
#' multiple SHAP analyses to identify universal vs. tier-specific features.
#'
#' @param shap_results List containing tier results from run_tiered_shap_analysis
#' @param n_top Integer. Number of top features to consider (default 25)
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{feature}{Feature name}
#'     \item{rank_<tier>}{Rank in each tier (NA if not in top n)}
#'     \item{mean_rank}{Mean rank across tiers where feature appears}
#'     \item{n_top_tiers}{Number of tiers where feature is in top n}
#'     \item{rank_variance}{Variance in ranks across tiers}
#'     \item{interpretation}{Classification of feature importance pattern}
#'   }
#'
#' @noRd
compare_tier_rankings <- function(shap_results, n_top = 25) {

  tier_names <- c("modern_quality", "modern_discord", "legacy_divergence", "robust_quality")
  available_tiers <- intersect(tier_names, names(shap_results))
  available_tiers <- available_tiers[vapply(available_tiers, function(t) {
    !is.null(shap_results[[t]])
  }, logical(1))]

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

  comparison$mean_rank <- rowMeans(comparison[, rank_cols, drop = FALSE], na.rm = TRUE)
  comparison$n_top_tiers <- rowSums(!is.na(comparison[, rank_cols, drop = FALSE]))
  comparison$rank_variance <- apply(comparison[, rank_cols, drop = FALSE], 1, var, na.rm = TRUE)

  # Sort by mean rank
  comparison <- comparison[order(comparison$mean_rank), ]

  # Add interpretation column
  comparison$interpretation <- vapply(1:nrow(comparison), function(i) {
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
  }, character(1))

  rownames(comparison) <- NULL
  return(comparison)
}


#' Print method for SHAP analysis results
#'
#' @param x A mutateR_shap_analysis object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.mutateR_shap_analysis <- function(x, ...) {

  cat("\n")
  cat("\u2554", paste(rep("\u2550", 63), collapse = ""), "\u2557\n", sep = "")
  cat("\u2551         mutateR SHAP Feature Importance Analysis              \u2551\n")
  cat("\u255A", paste(rep("\u2550", 63), collapse = ""), "\u255D\n\n", sep = "")

  cat("TARGET VARIABLE\n")
  cat(paste(rep("\u2500", 65), collapse = ""), "\n")
  cat("Target:          ", x$metadata$target, "\n")

  target_interp <- get_target_interpretation(x$metadata$target)
  if (!is.null(target_interp)) {
    cat("Interpretation:  ", target_interp, "\n")
  }
  cat("\n")

  cat("DATASET\n")
  cat(paste(rep("\u2500", 65), collapse = ""), "\n")
  cat("gRNAs analysed:  ", format(x$metadata$n_gRNAs, big.mark = ","), "\n")
  cat("Genes included:  ", x$metadata$n_genes, "\n")
  cat("Features:        ", x$metadata$n_features, "\n")
  cat("Modern methods:  ", paste(x$metadata$modern_methods, collapse = ", "), "\n")
  cat("Legacy methods:  ", paste(x$metadata$legacy_methods, collapse = ", "), "\n")
  cat("\n")

  cat("MODEL PERFORMANCE\n")
  cat(paste(rep("\u2500", 65), collapse = ""), "\n")
  cat(sprintf("R\u00B2:              %.3f\n", x$performance$r_squared))
  cat(sprintf("RMSE:            %.3f\n", x$performance$rmse))
  cat(sprintf("Spearman \u03C1:      %.3f\n", x$performance$spearman_cor))
  cat(sprintf("XGBoost rounds:  %d\n", x$performance$n_rounds))
  cat("\n")

  cat("TOP 15 FEATURES (by mean |SHAP|)\n")
  cat(paste(rep("\u2500", 65), collapse = ""), "\n")

  top15 <- head(x$feature_importance, 15)
  top15_display <- data.frame(
    rank = top15$rank,
    feature = top15$feature,
    mean_abs_shap = sprintf("%.4f", top15$mean_abs_shap),
    direction = top15$direction,
    stringsAsFactors = FALSE
  )
  print(top15_display, row.names = FALSE)

  cat("\n")
  cat(paste(rep("\u2550", 63), collapse = ""), "\n")
  cat("Access plots:           $plots$beeswarm, $plots$bar\n")
  cat("Access SHAP values:     $shap_values\n")
  cat("Access full importance: $feature_importance\n")
  cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")

  invisible(x)
}


#' Print method for tiered SHAP analysis
#'
#' @param x A mutateR_tiered_shap object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.mutateR_tiered_shap <- function(x, ...) {

  cat("\n")
  cat("\u2554", paste(rep("\u2550", 63), collapse = ""), "\u2557\n", sep = "")
  cat("\u2551           mutateR Tiered SHAP Analysis Results                \u2551\n")
  cat("\u255A", paste(rep("\u2550", 63), collapse = ""), "\u255D\n\n", sep = "")

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
      cat(sprintf("\u2500\u2500\u2500 %s \u2500\u2500\u2500\n", tier_info[[tier]]))
      cat(sprintf("    gRNAs: %s | Features: %d | R\u00B2: %.3f\n",
                  format(x[[tier]]$metadata$n_gRNAs, big.mark = ","),
                  x[[tier]]$metadata$n_features,
                  x[[tier]]$performance$r_squared))

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
    cat("\u2500\u2500\u2500 Cross-Tier Feature Comparison \u2500\u2500\u2500\n")

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

  cat(paste(rep("\u2550", 63), collapse = ""), "\n")
  cat("Access individual tier results: $modern_quality, $modern_discord, etc.\n")
  cat("Access feature comparison: $comparison\n")
  cat("Access plots: $<tier>$plots$beeswarm or $<tier>$plots$bar\n")
  cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")

  invisible(x)
}


#' Summary method for SHAP analysis
#'
#' @param object A mutateR_shap_analysis object
#' @param n_features Integer. Number of top features to show (default 20)
#' @param ... Additional arguments (ignored)
#'
#' @export
summary.mutateR_shap_analysis <- function(object, n_features = 20, ...) {

  cat("\nSHAP Analysis Summary: ", object$metadata$target, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Performance: R\u00B2 = ", round(object$performance$r_squared, 3),
      ", \u03C1 = ", round(object$performance$spearman_cor, 3), "\n\n")

  cat("Top ", n_features, " Features:\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")

  top_features <- head(object$feature_importance, n_features)

  for (i in seq_len(nrow(top_features))) {
    cat(sprintf("%2d. %-25s %s (%.4f)\n",
                top_features$rank[i],
                top_features$feature[i],
                ifelse(top_features$direction[i] == "positive", "\u2191", "\u2193"),
                top_features$mean_abs_shap[i]))
  }

  invisible(object)
}
