#' @title SHAP Visualization Functions
#' @description
#' Visualization functions for SHAP analysis results including beeswarm plots,
#' bar plots, dependence plots, and cross-tier comparisons. Supports rasterization
#' for large datasets and publication-quality output.
#'
#' @name SHAP_visualization
NULL


#' Generate SHAP visualization plots
#'
#' Creates beeswarm and bar plots for SHAP feature importance with enhanced
#' annotations including sample size, model performance, and target interpretation.
#' Supports rasterization for large datasets to improve rendering performance.
#'
#' @param shap_values A shapviz object containing SHAP values
#' @param target Character. Target variable name for labeling
#' @param n_gRNAs Integer. Number of gRNAs in the analysis
#' @param n_genes Integer. Number of genes in the analysis
#' @param r_squared Numeric. Model R-squared value
#' @param spearman_cor Numeric. Spearman correlation between predictions and actual
#' @param modern_methods Character vector. Modern scoring methods used
#' @param legacy_methods Character vector. Legacy scoring methods used
#' @param n_top_features Integer. Number of top features to display (default 25)
#' @param rasterize Logical. Rasterize beeswarm points (default NULL = auto)
#' @param raster_dpi Integer. DPI for rasterization (default 300)
#' @param subsample_beeswarm Integer or NULL. Number of points to subsample for beeswarm
#'
#' @return List of ggplot objects:
#'   \describe{
#'     \item{beeswarm}{Beeswarm plot showing SHAP value distributions}
#'     \item{bar}{Bar plot showing mean |SHAP| values}
#'     \item{waterfall_best}{Waterfall plot for highest-scoring gRNA}
#'     \item{waterfall_worst}{Waterfall plot for lowest-scoring gRNA}
#'   }
#'
#' @details
#' For datasets with >5,000 gRNAs, rasterization is automatically enabled to
#' improve plot rendering performance. For datasets with >20,000 gRNAs,
#' stratified subsampling is applied to the beeswarm plot while maintaining
#' the distribution of SHAP values.
#'
#' Plot annotations include:
#' \itemize{
#'   \item Target variable interpretation
#'   \item Sample size (n gRNAs from n genes)
#'   \item Model performance (R², Spearman ρ)
#'   \item Method configuration (modern/legacy)
#'   \item Subsample note if applicable
#' }
#'
#' @seealso \code{\link{save_shap_plots}} to export plots,
#'   \code{\link{plot_shap_dependence}} for feature dependence plots
#'
#' @noRd
generate_shap_plots <- function(shap_values,
                                target,
                                n_gRNAs,
                                n_genes,
                                r_squared,
                                spearman_cor,
                                modern_methods,
                                legacy_methods,
                                n_top_features = 25,
                                rasterize = NULL,
                                raster_dpi = 300,
                                subsample_beeswarm = NULL) {

  # ═══════════════════════════════════════════════════════════════
  # Determine Rasterization and Subsampling
  # ═══════════════════════════════════════════════════════════════

  # Auto-determine rasterization based on sample size
  if (is.null(rasterize)) {
    rasterize <- n_gRNAs > 5000
  }

  # Auto-determine subsampling for very large datasets
  if (is.null(subsample_beeswarm)) {
    if (n_gRNAs > 50000) {
      subsample_beeswarm <- 20000
    } else if (n_gRNAs > 20000) {
      subsample_beeswarm <- 15000
    }
  }

  # Check for ggrastr package if rasterizing
  has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
  if (rasterize && !has_ggrastr) {
    warning("Package 'ggrastr' not available. Install with: install.packages('ggrastr')\n",
            "Falling back to non-rasterized plots (may be slow).")
    rasterize <- FALSE
  }

  # ═══════════════════════════════════════════════════════════════
  # Build Plot Annotations
  # ═══════════════════════════════════════════════════════════════

  target_interp <- get_target_interpretation(target)
  subtitle <- paste0(
    if (!is.null(target_interp)) paste0(target_interp, "\n") else "",
    format(n_gRNAs, big.mark = ","), " gRNAs from ", n_genes, " genes | ",
    "R\u00B2 = ", sprintf("%.3f", r_squared), ", \u03C1 = ", sprintf("%.3f", spearman_cor)
  )

  caption <- paste0(
    "Modern methods: ", paste(modern_methods, collapse = ", "),
    if (length(legacy_methods) > 0) {
      paste0(" | Legacy: ", paste(legacy_methods, collapse = ", "))
    } else ""
  )

  # Add subsample note if applicable
  subsample_note <- NULL
  if (!is.null(subsample_beeswarm) && n_gRNAs > subsample_beeswarm) {
    subsample_note <- paste0("Beeswarm: ", format(subsample_beeswarm, big.mark = ","),
                             " points shown (stratified subsample)")
  }

  plots <- list()

  # ═══════════════════════════════════════════════════════════════
  # Prepare SHAP Values for Beeswarm (with optional subsampling)
  # ═══════════════════════════════════════════════════════════════

  shap_for_beeswarm <- shap_values
  if (!is.null(subsample_beeswarm) && n_gRNAs > subsample_beeswarm) {
    shap_for_beeswarm <- subsample_shap_values(shap_values, subsample_beeswarm)
  }

  # ═══════════════════════════════════════════════════════════════
  # Beeswarm Plot
  # ═══════════════════════════════════════════════════════════════

  plots$beeswarm <- tryCatch({

    # Create base beeswarm
    p_bee <- shapviz::sv_importance(
      shap_for_beeswarm,
      kind = "beeswarm",
      max_display = n_top_features,
      alpha = 0.5,
      size = 0.8
    )

    # Apply rasterization if requested and available
    # Note: ggrastr::rasterise() wraps the plot, not added with +
    if (rasterize && has_ggrastr) {
      p_bee <- ggrastr::rasterise(p_bee, dpi = raster_dpi)
    }

    # Build caption with subsample note
    beeswarm_caption <- caption
    if (!is.null(subsample_note)) {
      beeswarm_caption <- paste0(caption, "\n", subsample_note)
    }

    # Add enhanced labels
    p_bee <- p_bee +
      ggplot2::labs(
        title = paste0("SHAP Feature Importance: ", target),
        subtitle = subtitle,
        x = "SHAP value (impact on prediction)",
        y = NULL,
        caption = beeswarm_caption
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey30", hjust = 0.5,
                                              lineheight = 1.2),
        plot.caption = ggplot2::element_text(size = 8, color = "grey50", hjust = 0),
        legend.position = "right"
      )

    p_bee

  }, error = function(e) {
    warning("Beeswarm plot generation failed: ", e$message)
    NULL
  })

  # ═══════════════════════════════════════════════════════════════
  # Bar Plot
  # ═══════════════════════════════════════════════════════════════

  plots$bar <- tryCatch({

    shapviz::sv_importance(
      shap_values,
      kind = "bar",
      max_display = n_top_features,
      fill = "#3182bd"
    ) +
      ggplot2::labs(
        title = paste0("Mean |SHAP| Importance: ", target),
        subtitle = subtitle,
        x = "Mean |SHAP value|",
        y = NULL,
        caption = caption
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey30", hjust = 0.5,
                                              lineheight = 1.2),
        plot.caption = ggplot2::element_text(size = 8, color = "grey50", hjust = 0)
      )

  }, error = function(e) {
    warning("Bar plot generation failed: ", e$message)
    NULL
  })

  # ═══════════════════════════════════════════════════════════════
  # Waterfall Plots for Extreme Examples
  # ═══════════════════════════════════════════════════════════════

  plots$waterfall_best <- tryCatch({
    shap_matrix <- shapviz::get_shap_values(shap_values)
    baseline <- attr(shap_values, "baseline")
    if (is.null(baseline)) baseline <- 0
    predictions <- rowSums(shap_matrix) + baseline
    best_idx <- which.max(predictions)

    shapviz::sv_waterfall(shap_values, row_id = best_idx, max_display = 12) +
      ggplot2::labs(
        title = "SHAP Waterfall: Highest-Scoring gRNA",
        subtitle = "Feature contributions to prediction"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, color = "grey30", hjust = 0.5)
      )
  }, error = function(e) NULL)

  plots$waterfall_worst <- tryCatch({
    shap_matrix <- shapviz::get_shap_values(shap_values)
    baseline <- attr(shap_values, "baseline")
    if (is.null(baseline)) baseline <- 0
    predictions <- rowSums(shap_matrix) + baseline
    worst_idx <- which.min(predictions)

    shapviz::sv_waterfall(shap_values, row_id = worst_idx, max_display = 12) +
      ggplot2::labs(
        title = "SHAP Waterfall: Lowest-Scoring gRNA",
        subtitle = "Feature contributions to prediction"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, color = "grey30", hjust = 0.5)
      )
  }, error = function(e) NULL)

  return(plots)
}


#' Subsample SHAP values with stratification
#'
#' Performs stratified subsampling that maintains the distribution of SHAP values
#' for visualization of very large datasets. Stratification is based on overall
#' SHAP magnitude to ensure representation across the full range of predictions.
#'
#' @param shap_values A shapviz object
#' @param n_subsample Integer. Target number of samples
#' @param seed Integer. Random seed for reproducibility (default 42)
#'
#' @return A subsampled shapviz object
#'
#' @noRd
subsample_shap_values <- function(shap_values, n_subsample, seed = 42) {

  shap_matrix <- shapviz::get_shap_values(shap_values)
  X <- shapviz::get_feature_values(shap_values)
  n_total <- nrow(shap_matrix)

  if (n_total <= n_subsample) {
    return(shap_values)
  }

  set.seed(seed)

  # Stratified sampling based on overall SHAP magnitude
  total_shap <- rowSums(abs(shap_matrix))

  # Create quantile-based strata
  n_strata <- 10
  strata <- cut(
    total_shap,
    breaks = quantile(total_shap, probs = seq(0, 1, length.out = n_strata + 1)),
    labels = FALSE,
    include.lowest = TRUE
  )

  # Sample proportionally from each stratum
  samples_per_stratum <- ceiling(n_subsample / n_strata)

  selected_idx <- unlist(lapply(1:n_strata, function(s) {
    stratum_idx <- which(strata == s)
    n_select <- min(samples_per_stratum, length(stratum_idx))
    if (n_select > 0) {
      sample(stratum_idx, n_select)
    } else {
      integer(0)
    }
  }))

  # Trim to exact size if needed
  if (length(selected_idx) > n_subsample) {
    selected_idx <- sample(selected_idx, n_subsample)
  }

  selected_idx <- sort(selected_idx)

  # Subset the data
  X_sub <- X[selected_idx, , drop = FALSE]
  S_sub <- shap_matrix[selected_idx, , drop = FALSE]

  # Get baseline
  baseline <- attr(shap_values, "baseline")
  if (is.null(baseline)) baseline <- 0

  # Create new shapviz object
  shapviz::shapviz(S_sub, X = X_sub, baseline = baseline)
}


#' Save SHAP plots in multiple formats
#'
#' Convenience function to save plots in publication-ready formats.
#'
#' @param plots List of ggplot objects from SHAP analysis (e.g., result$plots)
#' @param output_dir Character. Output directory path
#' @param prefix Character. File name prefix (default "shap")
#' @param formats Character vector. Output formats (default: c("png", "pdf"))
#' @param width Numeric. Plot width in inches (default 10)
#' @param height Numeric. Plot height in inches (default 8)
#' @param dpi Integer. Resolution for raster formats (default 300)
#'
#' @return Character vector of saved file paths (invisibly)
#'
#' @examples
#' \dontrun{
#' results <- analyze_feature_importance(batch_result)
#'
#' # Save all plots
#' save_shap_plots(results$plots, "figures", prefix = "quality_modern")
#'
#' # Save only specific plots
#' save_shap_plots(
#'   results$plots[c("beeswarm", "bar")],
#'   "figures",
#'   prefix = "main_figure",
#'   formats = "pdf",
#'   width = 12,
#'   height = 10
#' )
#' }
#'
#' @export
save_shap_plots <- function(plots,
                            output_dir,
                            prefix = "shap",
                            formats = c("png", "pdf"),
                            width = 10,
                            height = 8,
                            dpi = 300) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  saved_files <- character()

  for (plot_name in names(plots)) {
    if (is.null(plots[[plot_name]])) next

    for (fmt in formats) {
      filename <- file.path(output_dir, paste0(prefix, "_", plot_name, ".", fmt))

      tryCatch({
        ggplot2::ggsave(
          filename = filename,
          plot = plots[[plot_name]],
          width = width,
          height = height,
          dpi = dpi
        )
        saved_files <- c(saved_files, filename)
      }, error = function(e) {
        warning("Failed to save ", filename, ": ", e$message)
      })
    }
  }

  message("Saved ", length(saved_files), " plot files to: ", output_dir)
  invisible(saved_files)
}


#' Generate SHAP dependence plots for specific features
#'
#' Creates dependence plots showing the relationship between a feature's
#' value and its SHAP contribution, optionally colored by an interaction feature.
#'
#' @param shap_result Output from analyze_feature_importance()
#' @param features Character vector of feature names to plot. Default NULL uses top 4 features.
#' @param color_by Character. Feature to use for coloring points.
#'        Default "auto" for automatic selection based on interaction strength.
#' @param ncol Integer. Number of columns in plot grid (default 2)
#' @param rasterize Logical. Rasterize points for large datasets (default TRUE if n > 5000)
#' @param raster_dpi Integer. DPI for rasterization (default 300)
#'
#' @return A ggplot object (combined with patchwork if available) or list of ggplot objects
#'
#' @export
plot_shap_dependence <- function(shap_result,
                                 features = NULL,
                                 color_by = "auto",
                                 ncol = 2,
                                 rasterize = NULL,
                                 raster_dpi = 300) {

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

  # Auto-determine rasterization
  n_points <- nrow(shap_result$data$X)
  if (is.null(rasterize)) {
    rasterize <- n_points > 5000
  }

  has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
  if (rasterize && !has_ggrastr) {
    warning("Package 'ggrastr' not available for rasterization.")
    rasterize <- FALSE
  }

  # Generate plots
  plots <- lapply(features, function(feat) {
    tryCatch({
      p <- shapviz::sv_dependence(shap_result$shap_values, v = feat, color_var = color_by) +
        ggplot2::ggtitle(feat) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
        )

      # Rasterize by wrapping the plot, not adding with +
      if (rasterize && has_ggrastr) {
        p <- ggrastr::rasterise(p, dpi = raster_dpi)
      }

      p
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
    combined <- patchwork::wrap_plots(plots, ncol = ncol) +
      patchwork::plot_annotation(
        title = paste0("SHAP Dependence Plots: ", shap_result$metadata$target),
        subtitle = paste0(format(n_points, big.mark = ","), " gRNAs from ",
                          shap_result$metadata$n_genes, " genes"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 10, color = "grey30", hjust = 0.5)
        )
      )
    return(combined)
  }

  return(plots)
}


#' Compare SHAP results between two analyses
#'
#' Generates comparison visualizations and statistics for two SHAP analyses.
#' Useful for comparing feature importance between targets (e.g., modern quality
#' vs. legacy divergence) or between different datasets.
#'
#' @param shap1 First mutateR_shap_analysis object
#' @param shap2 Second mutateR_shap_analysis object
#' @param label1 Character. Label for first analysis (default: uses target name)
#' @param label2 Character. Label for second analysis (default: uses target name)
#' @param n_features Integer. Number of top features to compare (default 20)
#'
#' @return List containing:
#'   \describe{
#'     \item{comparison}{Data frame with feature rankings from both analyses}
#'     \item{plot}{ggplot scatter plot comparing ranks}
#'     \item{summary}{List with n_shared_top and rank_correlation}
#'   }
#'
#' @examples
#' \dontrun{
#' shap_quality <- analyze_feature_importance(batch, target = "quality_percentile_modern")
#' shap_discord <- analyze_feature_importance(batch, target = "discordance_modern")
#'
#' comp <- compare_shap_analyses(shap_quality, shap_discord)
#'
#' # View comparison table
#' head(comp$comparison)
#'
#' # View scatter plot
#' print(comp$plot)
#'
#' # Check rank correlation
#' comp$summary$rank_correlation
#' }
#'
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
  comparison$category <- vapply(1:nrow(comparison), function(i) {
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
  }, character(1))

  # Sort by minimum rank
  comparison$min_rank <- pmin(comparison[[rank_col1]], comparison[[rank_col2]], na.rm = TRUE)
  comparison <- comparison[order(comparison$min_rank), ]
  rownames(comparison) <- NULL

  # Create visualization
  plot_df <- head(comparison, n_features * 2)
  plot_df <- plot_df[plot_df$category != "Less important in both", ]

  rank_plot <- NULL

  if (nrow(plot_df) > 0) {

    rank_plot <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data[[rank_col1]],
        y = .data[[rank_col2]],
        color = category
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_text(
        ggplot2::aes(label = feature),
        size = 2.5, hjust = -0.1, vjust = -0.1,
        check_overlap = TRUE, show.legend = FALSE
      ) +
      ggplot2::scale_color_brewer(palette = "Set2") +
      ggplot2::labs(
        title = paste("Feature Rank Comparison:", label1, "vs", label2),
        x = paste("Rank in", label1),
        y = paste("Rank in", label2),
        color = "Category"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom"
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
  }

  # Calculate correlation
  rank_cor <- cor(
    comparison[[rank_col1]],
    comparison[[rank_col2]],
    use = "complete.obs",
    method = "spearman"
  )

  return(list(
    comparison = comparison,
    plot = rank_plot,
    summary = list(
      n_shared_top = sum(comparison$category == "Important in both", na.rm = TRUE),
      rank_correlation = rank_cor,
      label1 = label1,
      label2 = label2
    )
  ))
}


#' Create summary visualization across all SHAP tiers
#'
#' Creates a heatmap visualization showing feature importance rankings
#' across multiple SHAP analysis tiers (quality, discord, legacy, robust).
#'
#' @param tiered_shap Output from run_tiered_shap_analysis()
#' @param n_features Integer. Number of top features to show (default 20)
#' @param cluster_features Logical. Order features by mean importance across tiers (default TRUE)
#'
#' @return A ggplot heatmap object
#'
#' @examples
#' \dontrun{
#' tiered <- run_tiered_shap_analysis(batch_result)
#' plot_tiered_comparison(tiered)
#'
#' # Show more features
#' plot_tiered_comparison(tiered, n_features = 30)
#' }
#'
#' @export
plot_tiered_comparison <- function(tiered_shap,
                                   n_features = 20,
                                   cluster_features = TRUE) {

  if (!inherits(tiered_shap, "mutateR_tiered_shap")) {
    stop("Input must be a mutateR_tiered_shap object from run_tiered_shap_analysis()")
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Install with: install.packages('tidyr')")
  }

  comparison <- tiered_shap$comparison

  if (is.null(comparison) || nrow(comparison) == 0) {
    stop("No comparison data available in tiered SHAP results")
  }

  # Get rank columns
  rank_cols <- grep("^rank_", names(comparison), value = TRUE)

  if (length(rank_cols) < 2) {
    stop("Need at least 2 tiers for comparison plot")
  }

  # Filter to top features (those that appear in top N of any tier)
  comparison$in_any_top <- apply(comparison[, rank_cols, drop = FALSE], 1, function(x) {
    any(x <= n_features, na.rm = TRUE)
  })

  plot_data <- comparison[comparison$in_any_top, ]

  if (nrow(plot_data) == 0) {
    warning("No features found in top ", n_features, " of any tier")
    return(NULL)
  }

  # Convert to long format for plotting
  plot_long <- tidyr::pivot_longer(
    plot_data,
    cols = tidyr::all_of(rank_cols),
    names_to = "tier",
    values_to = "rank",
    names_prefix = "rank_"
  )

  # Clean up tier names for display
  tier_labels <- c(
    "modern_quality" = "Quality\n(Modern)",
    "modern_discord" = "Discord\n(Modern)",
    "legacy_divergence" = "Legacy\nDivergence",
    "robust_quality" = "Quality\n(Robust)"
  )

  plot_long$tier_label <- tier_labels[plot_long$tier]
  plot_long$tier_label[is.na(plot_long$tier_label)] <- plot_long$tier[is.na(plot_long$tier_label)]

  # Order features
  if (cluster_features) {
    feature_order <- plot_data$feature[order(plot_data$mean_rank)]
  } else {
    first_tier <- rank_cols[1]
    feature_order <- plot_data$feature[order(plot_data[[first_tier]])]
  }

  plot_long$feature <- factor(plot_long$feature, levels = rev(feature_order))

  # Create heatmap
  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = tier_label, y = feature, fill = rank)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(rank), "", as.character(rank))),
      size = 3, color = "white"
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      direction = -1,
      na.value = "grey90",
      limits = c(1, n_features),
      oob = scales::squish,
      name = "Importance\nRank"
    ) +
    ggplot2::labs(
      title = "Feature Importance Across SHAP Tiers",
      subtitle = paste0("Top ", n_features, " features from each tier | Lower rank = more important"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30", hjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  return(p)
}


#' Create feature importance bar plot with category coloring
#'
#' Creates a bar plot of feature importance with bars colored by feature category
#' (compositional, positional, dinucleotide, etc.).
#'
#' @param shap_result Output from analyze_feature_importance()
#' @param n_features Integer. Number of top features to display (default 20)
#' @param show_direction Logical. Show arrows indicating effect direction (default TRUE)
#'
#' @return A ggplot object
#'
#' @export
plot_importance_by_category <- function(shap_result,
                                        n_features = 20,
                                        show_direction = TRUE) {

  if (!inherits(shap_result, "mutateR_shap_analysis")) {
    stop("Input must be a mutateR_shap_analysis object")
  }

  # Get feature categories
  categories <- get_feature_categories()

  # Create category lookup
  category_lookup <- character()
  for (cat_name in names(categories)) {
    for (feat in categories[[cat_name]]) {
      category_lookup[feat] <- cat_name
    }
  }

  # Get top features
  plot_data <- head(shap_result$feature_importance, n_features)

  # Add category
  plot_data$category <- vapply(plot_data$feature, function(f) {
    if (f %in% names(category_lookup)) {
      return(category_lookup[f])
    }
    return("other")
  }, character(1))

  # Add direction symbol BEFORE creating the plot
  if (show_direction) {
    plot_data$direction_symbol <- ifelse(
      plot_data$direction == "positive",
      "\u2191",  # up arrow
      "\u2193"   # down arrow
    )
  }

  # Order features by importance
  plot_data$feature <- factor(plot_data$feature, levels = rev(plot_data$feature))

  # Define category colors
  category_colors <- c(
    "composition" = "#1b9e77",
    "homopolymer" = "#d95f02",
    "dinucleotide" = "#7570b3",
    "positional" = "#e7298a",
    "thermodynamic" = "#66a61e",
    "structural" = "#e6ab02",
    "context" = "#a6761d",
    "other" = "#666666"
  )

  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = mean_abs_shap, y = feature, fill = category)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::scale_fill_manual(values = category_colors, name = "Feature Category") +
    ggplot2::labs(
      title = paste0("Feature Importance by Category: ", shap_result$metadata$target),
      subtitle = paste0(format(shap_result$metadata$n_gRNAs, big.mark = ","),
                        " gRNAs | R\u00B2 = ", round(shap_result$performance$r_squared, 3)),
      x = "Mean |SHAP value|",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30", hjust = 0.5),
      legend.position = "right"
    )

  # Add direction arrows if requested
  if (show_direction) {
    p <- p +
      ggplot2::geom_text(
        data = plot_data,  # Explicitly pass data with direction_symbol
        ggplot2::aes(x = mean_abs_shap, y = feature, label = direction_symbol),
        hjust = -0.5,
        size = 4,
        inherit.aes = FALSE  # Don't inherit fill aesthetic
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))
  }

  return(p)
}


#' Create waterfall plot for a specific gRNA
#'
#' Creates a SHAP waterfall plot showing feature contributions for a specific gRNA,
#' identified by index or sequence.
#'
#' @param shap_result Output from analyze_feature_importance()
#' @param index Integer. Row index of gRNA to plot (default NULL)
#' @param sequence Character. Protospacer sequence to match (default NULL)
#' @param max_display Integer. Maximum features to show (default 12)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' results <- analyze_feature_importance(batch_result)
#'
#' # Plot by index
#' plot_waterfall_for_grna(results, index = 1)
#'
#' # Plot by sequence
#' plot_waterfall_for_grna(results, sequence = "AGCTGATCGATCGATCGATC")
#' }
#'
#' @export
plot_waterfall_for_grna <- function(shap_result,
                                    index = NULL,
                                    sequence = NULL,
                                    max_display = 12) {

  if (!inherits(shap_result, "mutateR_shap_analysis")) {
    stop("Input must be a mutateR_shap_analysis object")
  }

  if (is.null(index) && is.null(sequence)) {
    stop("Must provide either 'index' or 'sequence'")
  }

  # Find index by sequence if needed
  if (!is.null(sequence)) {
    if (!"protospacer_sequence" %in% names(shap_result$data$features)) {
      stop("Cannot match by sequence: protospacer_sequence not in features data")
    }

    matches <- which(shap_result$data$features$protospacer_sequence == toupper(sequence))

    if (length(matches) == 0) {
      stop("Sequence not found in data")
    }
    if (length(matches) > 1) {
      warning("Multiple matches found, using first")
    }

    index <- matches[1]
  }

  # Validate index
  n_rows <- nrow(shap_result$data$X)
  if (index < 1 || index > n_rows) {
    stop("Index out of range. Valid range: 1 to ", n_rows)
  }

  # Create waterfall plot
  p <- shapviz::sv_waterfall(shap_result$shap_values, row_id = index, max_display = max_display) +
    ggplot2::labs(
      title = paste0("SHAP Waterfall: gRNA #", index),
      subtitle = paste0("Target: ", shap_result$metadata$target)
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey30", hjust = 0.5)
    )

  return(p)
}
