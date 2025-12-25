#' @title SHAP Ablation Analysis
#' @description
#' Ablation analysis functions for understanding gRNA features in different
#' experimental contexts. Enables removal of transcription-related features
#' to simulate synthetic gRNA / RNP delivery contexts where Pol III
#' transcription is not a factor.
#'
#' @name SHAP_ablation
NULL


#' Get list of transcription-related features
#'
#' Returns feature names that are related to gRNA transcription efficiency
#' rather than intrinsic Cas9:gRNA:DNA recognition. These features are
#' relevant for plasmid/lentiviral delivery but not for synthetic gRNA/RNP delivery.
#'
#' @param include_extended Logical. Include extended 5' region features (positions 2-3).
#'        Default TRUE.
#'
#' @return Character vector of transcription-related feature names
#'
#' @details
#' Features are classified as transcription-related if they:
#' \itemize{
#'   \item Affect Pol III termination (TT dinucleotides, poly-T runs)
#'   \item Affect transcription initiation (position 1 nucleotide for U6/T7)
#'   \item Affect transcript processing/stability (5' region)
#' }
#'
#' \strong{Core transcription features:}
#' \itemize{
#'   \item \code{freq_TT}: TT dinucleotide frequency (Pol III terminator signal)
#'   \item \code{has_polyT}: Presence of TTTT run
#'   \item \code{pos1_*}: Position 1 nucleotides (U6 +1 preference for G)
#'   \item \code{max_homopolymer}: Longest homopolymer run
#' }
#'
#' \strong{Extended features (if include_extended = TRUE):}
#' \itemize{
#'   \item \code{pos2_*}, \code{pos3_*}: Extended 5' region affecting processing
#'   \item \code{has_polyA/G/C}: Other homopolymer runs
#' }
#'
#' @examples
#' # Get core transcription features
#' tx_core <- get_transcription_features(include_extended = FALSE)
#' length(tx_core)
#'
#' # Get all transcription features
#' tx_all <- get_transcription_features()
#' length(tx_all)
#'
#' # Use for feature filtering
#' all_features <- get_feature_names()
#' intrinsic_features <- setdiff(all_features, tx_all)
#'
#' @seealso
#' \code{\link{run_transcription_ablation}} for ablation analysis,
#' \code{\link{get_feature_categories}} for all feature categories
#'
#' @export
get_transcription_features <- function(include_extended = TRUE) {

  # Core transcription features (always included)
  core_features <- c(
    # Pol III termination signals
    "freq_TT",
    "has_polyT",

    # Position 1 (U6/T7 +1 preference)
    "pos1_A",
    "pos1_T",
    "pos1_G",
    "pos1_C",

    # Homopolymer runs (affect transcription)
    "max_homopolymer"
  )

  if (!include_extended) {
    return(core_features)
  }

  # Extended features (5' region and other homopolymers)
  extended_features <- c(
    # Other homopolymer indicators
    "has_polyA",
    "has_polyG",
    "has_polyC",

    # 5' region (positions 2-3 affect processing)
    "pos2_A",
    "pos2_T",
    "pos2_G",
    "pos2_C",
    "pos3_A",
    "pos3_T",
    "pos3_G",
    "pos3_C"
  )

  return(c(core_features, extended_features))
}


#' Get list of intrinsic recognition features
#'
#' Returns feature names related to intrinsic Cas9:gRNA:DNA recognition,
#' excluding transcription-related features. These features should be
#' predictive regardless of gRNA delivery method.
#'
#' @return Character vector of intrinsic recognition feature names
#'
#' @details
#' Intrinsic recognition features include:
#' \itemize{
#'   \item Seed region nucleotides (positions 12-20)
#'   \item PAM-proximal preferences (positions 17-20)
#'   \item GC content (overall and regional)
#'   \item Thermodynamic properties
#'   \item Self-complementarity
#' }
#'
#' @examples
#' intrinsic <- get_intrinsic_features()
#' length(intrinsic)
#'
#' @seealso \code{\link{get_transcription_features}} for transcription features
#'
#' @export
get_intrinsic_features <- function() {

  all_features <- get_feature_names()
  tx_features <- get_transcription_features(include_extended = TRUE)

  setdiff(all_features, tx_features)
}


#' Run ablation analysis removing transcription features
#'
#' Performs SHAP analysis with transcription-related features removed,
#' simulating the feature importance profile for synthetic gRNA delivery.
#' This helps identify which features matter for intrinsic Cas9:gRNA:DNA
#' recognition, separate from transcription efficiency.
#'
#' @param batch_result A mutateR_consensus_batch object
#' @param target Character. Target variable to model.
#'        Default "quality_percentile_modern".
#' @param modern_methods Character vector of modern scoring methods.
#'        Default: c("deepspcas9", "deephf", "ruleset3")
#' @param legacy_methods Character vector of legacy scoring methods.
#'        Default: c("ruleset1")
#' @param n_top_features Integer. Number of top features to display (default 20)
#' @param include_extended Logical. Remove extended 5' features too (default TRUE)
#' @param rasterize Logical. Rasterize beeswarm plots (default NULL = auto)
#' @param raster_dpi Integer. DPI for rasterization (default 300)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return A mutateR_shap_analysis object with additional ablation metadata:
#'   \describe{
#'     \item{metadata$ablation$type}{"transcription"}
#'     \item{metadata$ablation$features_removed}{Character vector of removed features}
#'     \item{metadata$ablation$context}{Description of simulated context}
#'   }
#'
#' @details
#' The ablation analysis addresses a key limitation of existing training data:
#' most large-scale CRISPR screens use plasmid or lentiviral delivery where
#' gRNA expression from U6/H1 promoters is a major determinant of activity.
#'
#' Features like \code{freq_TT} (Pol III termination), \code{pos1_C} (poor U6
#' initiation), and \code{max_homopolymer} dominate standard SHAP analyses
#' because they affect transcription, not intrinsic Cas9 activity.
#'
#' By removing these features, we can ask: "What would predict gRNA efficacy
#' if transcription weren't a factor?" This is directly relevant for:
#' \itemize{
#'   \item Synthetic gRNA delivery
#'   \item RNP (ribonucleoprotein) delivery
#'   \item In vitro Cas9 activity assays
#'   \item Therapeutic applications avoiding DNA delivery
#' }
#'
#' @section Comparison with Full Analysis:
#' After running ablation, compare with full analysis:
#' \preformatted{
#' full_shap <- analyze_feature_importance(batch, target = "quality_percentile_modern")
#' ablated_shap <- run_transcription_ablation(batch, target = "quality_percentile_modern")
#'
#' comp <- compare_shap_analyses(full_shap, ablated_shap,
#'                               label1 = "Full", label2 = "Intrinsic")
#' print(comp$plot)
#' }
#'
#' @section Expected Changes:
#' Removing transcription features typically:
#' \itemize{
#'   \item Reduces model R² (transcription features explain substantial variance)
#'   \item Elevates seed region features (pos14-20) in importance
#'   \item Increases relative importance of GC content features
#'   \item May reveal thermodynamic/structural features previously masked
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
#' # Full analysis
#' full_shap <- analyze_feature_importance(batch)
#'
#' # Ablation analysis (synthetic gRNA context)
#' ablated_shap <- run_transcription_ablation(batch)
#'
#' # Compare
#' print(full_shap)
#' print(ablated_shap)
#'
#' # Detailed comparison
#' comp <- compare_shap_analyses(full_shap, ablated_shap,
#'                               label1 = "Plasmid", label2 = "Synthetic")
#' comp$plot
#' }
#'
#' @seealso
#' \code{\link{get_transcription_features}} for feature list,
#' \code{\link{compare_shap_analyses}} for comparing results,
#' \code{\link{run_custom_ablation}} for custom feature ablation
#'
#' @export
run_transcription_ablation <- function(batch_result,
                                       target = "quality_percentile_modern",
                                       modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                       legacy_methods = c("ruleset1"),
                                       n_top_features = 20,
                                       include_extended = TRUE,
                                       rasterize = NULL,
                                       raster_dpi = 300,
                                       quiet = FALSE) {

  if (!quiet) {
    cat("\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")
    cat("  Transcription Feature Ablation Analysis\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")
    cat("\nRemoving transcription-related features to simulate\n")
    cat("synthetic gRNA / RNP delivery context.\n\n")
  }

  tx_features <- get_transcription_features(include_extended = include_extended)

  if (!quiet) {
    cat("Features being removed (", length(tx_features), "):\n", sep = "")

    # Group for display
    core <- tx_features[tx_features %in% c("freq_TT", "has_polyT", "max_homopolymer")]
    pos1 <- tx_features[grepl("^pos1_", tx_features)]
    extended <- setdiff(tx_features, c(core, pos1))

    if (length(core) > 0) {
      cat("  Core: ", paste(core, collapse = ", "), "\n", sep = "")
    }
    if (length(pos1) > 0) {
      cat("  Position 1: ", paste(pos1, collapse = ", "), "\n", sep = "")
    }
    if (length(extended) > 0) {
      cat("  Extended: ", paste(extended, collapse = ", "), "\n", sep = "")
    }
    cat("\n")
  }

  # Run ablated analysis
  result <- analyze_feature_importance_ablated(
    batch_result = batch_result,
    target = target,
    exclude_features = tx_features,
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    rasterize = rasterize,
    raster_dpi = raster_dpi,
    quiet = quiet
  )

  # Add ablation metadata
  result$metadata$ablation <- list(
    type = "transcription",
    features_removed = tx_features,
    n_features_removed = length(tx_features),
    include_extended = include_extended,
    context = "Simulates synthetic gRNA / RNP delivery (no Pol III transcription)"
  )

  if (!quiet) {
    cat("\nAblation analysis complete.\n")
    cat("Original features: ", result$metadata$n_features + length(tx_features), "\n", sep = "")
    cat("Features removed:  ", length(tx_features), "\n", sep = "")
    cat("Features used:     ", result$metadata$n_features, "\n", sep = "")
    cat("\nCompare with full analysis using compare_shap_analyses()\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")
  }

  return(result)
}


#' Run custom ablation analysis
#'
#' Performs SHAP analysis with a custom set of features removed.
#' Useful for testing specific hypotheses about feature importance.
#'
#' @param batch_result A mutateR_consensus_batch object
#' @param exclude_features Character vector. Features to exclude from analysis.
#' @param target Character. Target variable to model.
#' @param modern_methods Character vector of modern scoring methods.
#' @param legacy_methods Character vector of legacy scoring methods.
#' @param n_top_features Integer. Number of top features to display (default 20)
#' @param ablation_name Character. Name for this ablation (for metadata).
#' @param ablation_description Character. Description of ablation rationale.
#' @param rasterize Logical. Rasterize beeswarm plots (default NULL = auto)
#' @param raster_dpi Integer. DPI for rasterization (default 300)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return A mutateR_shap_analysis object with ablation metadata
#'
#' @examples
#' \dontrun{
#' # Remove all positional features to test compositional-only model
#' positional <- get_feature_categories()$positional
#' ablated <- run_custom_ablation(
#'   batch,
#'   exclude_features = positional,
#'   ablation_name = "no_positional",
#'   ablation_description = "Compositional features only"
#' )
#'
#' # Remove seed region features
#' seed_features <- paste0("pos", 12:20, "_", rep(c("A","T","G","C"), each = 9))
#' ablated_seed <- run_custom_ablation(
#'   batch,
#'   exclude_features = seed_features,
#'   ablation_name = "no_seed",
#'   ablation_description = "Without seed region positional features"
#' )
#' }
#'
#' @export
run_custom_ablation <- function(batch_result,
                                exclude_features,
                                target = "quality_percentile_modern",
                                modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                legacy_methods = c("ruleset1"),
                                n_top_features = 20,
                                ablation_name = "custom",
                                ablation_description = NULL,
                                rasterize = NULL,
                                raster_dpi = 300,
                                quiet = FALSE) {

  if (length(exclude_features) == 0) {
    stop("exclude_features cannot be empty")
  }

  if (!quiet) {
    cat("\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")
    cat("  Custom Feature Ablation: ", ablation_name, "\n", sep = "")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")

    if (!is.null(ablation_description)) {
      cat("\n", ablation_description, "\n", sep = "")
    }

    cat("\nRemoving ", length(exclude_features), " features:\n", sep = "")
    if (length(exclude_features) <= 10) {
      cat("  ", paste(exclude_features, collapse = ", "), "\n", sep = "")
    } else {
      cat("  ", paste(head(exclude_features, 10), collapse = ", "),
          ", ... (", length(exclude_features) - 10, " more)\n", sep = "")
    }
    cat("\n")
  }

  # Run ablated analysis
  result <- analyze_feature_importance_ablated(
    batch_result = batch_result,
    target = target,
    exclude_features = exclude_features,
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    rasterize = rasterize,
    raster_dpi = raster_dpi,
    quiet = quiet
  )

  # Add ablation metadata
  result$metadata$ablation <- list(
    type = "custom",
    name = ablation_name,
    description = ablation_description,
    features_removed = exclude_features,
    n_features_removed = length(exclude_features)
  )

  if (!quiet) {
    cat("\nAblation analysis complete.\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")
  }

  return(result)
}


#' Internal function for ablated SHAP analysis
#'
#' Core implementation of feature ablation for SHAP analysis.
#' Called by run_transcription_ablation() and run_custom_ablation().
#'
#' @param batch_result A mutateR_consensus_batch object
#' @param target Target variable name
#' @param exclude_features Features to exclude
#' @param modern_methods Modern method names
#' @param legacy_methods Legacy method names
#' @param n_top_features Number of top features for plots
#' @param rasterize Whether to rasterize plots
#' @param raster_dpi DPI for rasterization
#' @param quiet Suppress messages
#'
#' @return A mutateR_shap_analysis object
#'
#' @noRd
analyze_feature_importance_ablated <- function(batch_result,
                                               target,
                                               exclude_features,
                                               modern_methods,
                                               legacy_methods,
                                               n_top_features,
                                               rasterize = NULL,
                                               raster_dpi = 300,
                                               quiet = FALSE) {

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required.")
  }

  if (!requireNamespace("shapviz", quietly = TRUE)) {
    stop("Package 'shapviz' is required.")
  }

  # ═══════════════════════════════════════════════════════════════
  # Handle Single Analysis Object
  # ═══════════════════════════════════════════════════════════════

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

  # ═══════════════════════════════════════════════════════════════
  # Aggregate Data
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Aggregating data across genes...")

  all_features <- list()
  all_targets <- list()

  for (gene in names(batch_result)) {
    obj <- batch_result[[gene]]
    if (is.null(obj$scores) || nrow(obj$scores) == 0) next

    feat <- extract_grna_features(obj$scores)
    feat$gene <- gene

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

  # ═══════════════════════════════════════════════════════════════
  # Prepare Data WITH ABLATION
  # ═══════════════════════════════════════════════════════════════

  feature_cols <- setdiff(names(features_df), c("grna_id", "gene"))

  # Check which excluded features actually exist
  existing_excluded <- intersect(exclude_features, feature_cols)
  missing_excluded <- setdiff(exclude_features, feature_cols)

  if (length(missing_excluded) > 0 && !quiet) {
    message("Note: ", length(missing_excluded), " excluded features not found in data")
  }

  # Remove excluded features
  feature_cols <- setdiff(feature_cols, existing_excluded)

  if (length(feature_cols) == 0) {
    stop("No features remaining after ablation")
  }

  if (!quiet) {
    message("Using ", length(feature_cols), " features (",
            length(existing_excluded), " removed)")
  }

  X <- as.matrix(features_df[, feature_cols])
  y <- targets_df[[target]]

  # Remove rows with NA
  valid <- complete.cases(X) & !is.na(y)
  X <- X[valid, , drop = FALSE]
  y <- y[valid]
  features_df_valid <- features_df[valid, ]
  targets_df_valid <- targets_df[valid, ]

  n_gRNAs <- nrow(X)

  if (!quiet) {
    message("Training on ", format(n_gRNAs, big.mark = ","), " observations")
  }

  # ═══════════════════════════════════════════════════════════════
  # Train XGBoost Model
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Training gradient boosting model...")

  dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

  params <- list(
    objective = "reg:squarederror",
    max_depth = 6,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 5
  )

  # Cross-validation for early stopping
  cv_result <- tryCatch({
    xgboost::xgb.cv(
      params = params,
      data = dtrain,
      nrounds = 500,
      nfold = 5,
      early_stopping_rounds = 20,
      verbose = 0
    )
  }, error = function(e) NULL)

  best_nrounds <- if (!is.null(cv_result)) cv_result$best_iteration else 100
  if (is.null(best_nrounds) || best_nrounds < 10) best_nrounds <- 100

  model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )

  # ═══════════════════════════════════════════════════════════════
  # Compute SHAP Values
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Computing SHAP values...")

  shap_values <- shapviz::shapviz(model, X_pred = X, X = X)

  # ═══════════════════════════════════════════════════════════════
  # Model Performance
  # ═══════════════════════════════════════════════════════════════

  predictions <- predict(model, X)

  ss_res <- sum((y - predictions)^2)
  ss_tot <- sum((y - mean(y))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((y - predictions)^2))
  pred_cor <- cor(y, predictions, method = "spearman")

  # ═══════════════════════════════════════════════════════════════
  # Generate Plots
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) message("Generating visualizations...")

  plots <- generate_shap_plots(
    shap_values = shap_values,
    target = paste0(target, " (Ablated)"),
    n_gRNAs = n_gRNAs,
    n_genes = n_genes,
    r_squared = r_squared,
    spearman_cor = pred_cor,
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    rasterize = rasterize,
    raster_dpi = raster_dpi
  )

  # ═══════════════════════════════════════════════════════════════
  # Feature Importance
  # ═══════════════════════════════════════════════════════════════

  shap_matrix <- shapviz::get_shap_values(shap_values)
  mean_abs_shap <- colMeans(abs(shap_matrix))

  importance_df <- data.frame(
    feature = names(mean_abs_shap),
    mean_abs_shap = as.numeric(mean_abs_shap),
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(-importance_df$mean_abs_shap), ]
  importance_df$rank <- seq_len(nrow(importance_df))

  mean_shap <- colMeans(shap_matrix)
  importance_df$mean_shap <- mean_shap[importance_df$feature]
  importance_df$direction <- ifelse(importance_df$mean_shap > 0, "positive", "negative")
  rownames(importance_df) <- NULL

  # ═══════════════════════════════════════════════════════════════
  # Return Result
  # ═══════════════════════════════════════════════════════════════

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
      excluded_features = existing_excluded
    )
  )

  class(result) <- c("mutateR_shap_analysis", class(result))

  if (!quiet) {
    message("Ablated model R\u00B2 = ", round(r_squared, 3),
            ", Spearman \u03C1 = ", round(pred_cor, 3))
  }

  return(result)
}


#' Compare full and ablated SHAP analyses
#'
#' Convenience function to run both full and ablated analyses and generate
#' a comprehensive comparison.
#'
#' @param batch_result A mutateR_consensus_batch object
#' @param target Character. Target variable to model.
#' @param ablation_type Character. Type of ablation: "transcription" (default),
#'        "positional", "dinucleotide", or "compositional".
#' @param modern_methods Character vector of modern scoring methods.
#' @param legacy_methods Character vector of legacy scoring methods.
#' @param n_top_features Integer. Number of top features to compare (default 20)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return List containing:
#'   \describe{
#'     \item{full}{Full SHAP analysis result}
#'     \item{ablated}{Ablated SHAP analysis result}
#'     \item{comparison}{Output from compare_shap_analyses()}
#'     \item{r_squared_change}{Change in R² from ablation}
#'     \item{summary}{Text summary of key findings}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- compare_full_vs_ablated(batch_result)
#'
#' # View R² change
#' result$r_squared_change
#'
#' # View comparison plot
#' result$comparison$plot
#'
#' # Print summary
#' cat(result$summary)
#' }
#'
#' @export
compare_full_vs_ablated <- function(batch_result,
                                    target = "quality_percentile_modern",
                                    ablation_type = c("transcription", "positional",
                                                      "dinucleotide", "compositional"),
                                    modern_methods = c("deepspcas9", "deephf", "ruleset3"),
                                    legacy_methods = c("ruleset1"),
                                    n_top_features = 20,
                                    quiet = FALSE) {

  ablation_type <- match.arg(ablation_type)

  if (!quiet) {
    cat("\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")
    cat("  Full vs Ablated SHAP Comparison\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Run Full Analysis
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) cat("Running full analysis...\n")

  full_shap <- analyze_feature_importance(
    batch_result,
    target = target,
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    quiet = TRUE
  )

  if (!quiet) {
    cat("  Full model R\u00B2: ", round(full_shap$performance$r_squared, 3), "\n\n", sep = "")
  }

  # ═══════════════════════════════════════════════════════════════
  # Run Ablated Analysis
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) cat("Running ablated analysis (", ablation_type, ")...\n", sep = "")

  # Determine features to exclude based on ablation type
  exclude_features <- switch(
    ablation_type,
    "transcription" = get_transcription_features(),
    "positional" = get_feature_categories()$positional,
    "dinucleotide" = get_feature_categories()$dinucleotide,
    "compositional" = get_feature_categories()$composition
  )

  ablation_description <- switch(
    ablation_type,
    "transcription" = "Simulates synthetic gRNA / RNP delivery",
    "positional" = "Tests importance of position-specific nucleotides",
    "dinucleotide" = "Tests importance of dinucleotide context",
    "compositional" = "Tests importance of GC content features"
  )

  ablated_shap <- run_custom_ablation(
    batch_result,
    exclude_features = exclude_features,
    target = target,
    modern_methods = modern_methods,
    legacy_methods = legacy_methods,
    n_top_features = n_top_features,
    ablation_name = ablation_type,
    ablation_description = ablation_description,
    quiet = TRUE
  )

  if (!quiet) {
    cat("  Ablated model R\u00B2: ", round(ablated_shap$performance$r_squared, 3), "\n\n", sep = "")
  }

  # ═══════════════════════════════════════════════════════════════
  # Compare Results
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) cat("Comparing results...\n")

  comparison <- compare_shap_analyses(
    full_shap,
    ablated_shap,
    label1 = "Full",
    label2 = "Ablated",
    n_features = n_top_features
  )

  # Calculate R² change
  r2_full <- full_shap$performance$r_squared
  r2_ablated <- ablated_shap$performance$r_squared
  r2_change <- r2_ablated - r2_full
  r2_pct_change <- 100 * r2_change / r2_full

  # ═══════════════════════════════════════════════════════════════
  # Generate Summary
  # ═══════════════════════════════════════════════════════════════

  # Find features that gained importance after ablation
  comp_df <- comparison$comparison
  comp_df$rank_change <- comp_df$rank_Full - comp_df$rank_Ablated
  gainers <- comp_df[!is.na(comp_df$rank_change) & comp_df$rank_change > 5, ]
  gainers <- gainers[order(-gainers$rank_change), ]

  # Find features unique to ablated top
  unique_ablated <- comp_df$feature[grepl("Unique to Ablated", comp_df$category)]

  summary_text <- sprintf(
    "Full vs Ablated Comparison Summary (%s ablation)
=============================================

Model Performance:
  Full model R²:     %.3f
  Ablated model R²:  %.3f
  Change:            %+.3f (%.1f%%)

Interpretation:
  Removing %d %s features %s the model's predictive power by %.1f%%.
  %s

Feature Shifts:
  Features shared in top %d: %d
  Rank correlation: %.3f

%s%s",
    ablation_type,
    r2_full,
    r2_ablated,
    r2_change,
    r2_pct_change,
    length(exclude_features),
    ablation_type,
    ifelse(r2_change < 0, "reduces", "increases"),
    abs(r2_pct_change),
    ifelse(r2_change < -0.05,
           "This suggests these features explain substantial variance in the scoring methods.",
           ifelse(r2_change > -0.02,
                  "This suggests these features have limited predictive importance.",
                  "This suggests moderate importance of these features.")),
    n_top_features,
    comparison$summary$n_shared_top,
    comparison$summary$rank_correlation,
    ifelse(nrow(gainers) > 0,
           paste0("\nFeatures gaining importance after ablation:\n  ",
                  paste(head(gainers$feature, 5), collapse = ", "), "\n"),
           ""),
    ifelse(length(unique_ablated) > 0,
           paste0("\nFeatures newly important after ablation:\n  ",
                  paste(head(unique_ablated, 5), collapse = ", "), "\n"),
           "")
  )

  if (!quiet) {
    cat("\n")
    cat(paste(rep("\u2550", 63), collapse = ""), "\n")
    cat(summary_text)
    cat(paste(rep("\u2550", 63), collapse = ""), "\n\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Return Results
  # ═══════════════════════════════════════════════════════════════

  return(list(
    full = full_shap,
    ablated = ablated_shap,
    comparison = comparison,
    r_squared_change = r2_change,
    r_squared_pct_change = r2_pct_change,
    ablation_type = ablation_type,
    excluded_features = exclude_features,
    summary = summary_text
  ))
}
