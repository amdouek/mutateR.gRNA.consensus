#' Visualise correlations between multiple scoring methods
#'
#' Generates a correlation matrix heatmap and optional pairwise scatter plots
#' to explore agreement between on-target scoring methods.
#'
#' @param multi_scores Data frame returned by \code{score_grnas_multi()}.
#' @param methods Character vector. Subset of methods to plot (default: all).
#' @param plot_type Character. One of:
#'   - "matrix": Correlation matrix heatmap (default)
#'   - "scatter": Pairwise scatter plots
#'   - "both": Both visualisations
#' @param cor_method Character. Correlation method: "spearman" (default) or "pearson".
#' @param highlight_top Integer. Number of top-ranked gRNAs to highlight in scatter plots (default 20).
#'
#' @return A ggplot object (or list of ggplot objects if plot_type = "both").
#'
#' @export
plot_score_correlations <- function(multi_scores,
                                    methods = NULL,
                                    plot_type = c("matrix", "scatter", "both"),
                                    cor_method = c("spearman", "pearson"),
                                    highlight_top = 20) {

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })

  plot_type <- match.arg(plot_type)
  cor_method <- match.arg(cor_method)

  # --- Identify score columns ---
  stored_methods <- attr(multi_scores, "methods")
  if (is.null(stored_methods)) {
    # Fallback: detect numeric columns that aren't metadata
    meta_cols <- c("grna_id", "protospacer_sequence", "sequence_context",
                   "exon_rank", "chr", "start", "end", "strand")
    stored_methods <- setdiff(names(multi_scores), meta_cols)
    stored_methods <- stored_methods[sapply(multi_scores[stored_methods], is.numeric)]
  }

  if (is.null(methods)) {
    methods <- stored_methods
  } else {
    invalid <- setdiff(methods, stored_methods)
    if (length(invalid) > 0) {
      stop("Methods not found in data: ", paste(invalid, collapse = ", "))
    }
  }

  if (length(methods) < 2) {
    stop("Need at least 2 methods to compute correlations.")
  }

  # --- Extract score matrix ---
  score_mat <- multi_scores[, methods, drop = FALSE]

  # --- Compute correlation matrix ---
  cor_mat <- cor(score_mat, use = "pairwise.complete.obs", method = cor_method)

  # --- Helper: Create correlation matrix heatmap ---
  plot_matrix <- function() {
    # Convert to long format for ggplot
    cor_df <- as.data.frame(as.table(cor_mat))
    names(cor_df) <- c("Method1", "Method2", "Correlation")

    # Create ordered factor for consistent display
    cor_df$Method1 <- factor(cor_df$Method1, levels = methods)
    cor_df$Method2 <- factor(cor_df$Method2, levels = rev(methods))

    p <- ggplot(cor_df, aes(x = Method1, y = Method2, fill = Correlation)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", Correlation)),
                color = ifelse(abs(cor_df$Correlation) > 0.5, "white", "black"),
                size = 4) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-1, 1),
                           name = paste0(tools::toTitleCase(cor_method), "\nCorrelation")) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      coord_fixed() +
      ggtitle(paste0("On-Target Score Correlations (", cor_method, ")"))

    return(p)
  }

  # --- Helper: Create pairwise scatter plots ---
  plot_scatter <- function() {
    # Generate all pairs
    pairs <- combn(methods, 2, simplify = FALSE)

    # Calculate ranks for highlighting
    ranks <- lapply(methods, function(m) {
      r <- rank(-score_mat[[m]], na.last = "keep", ties.method = "first")
      r
    })
    names(ranks) <- methods

    plot_list <- lapply(pairs, function(pair) {
      m1 <- pair[1]
      m2 <- pair[2]

      plot_df <- data.frame(
        x = score_mat[[m1]],
        y = score_mat[[m2]],
        rank_x = ranks[[m1]],
        rank_y = ranks[[m2]]
      )

      # Identify top candidates in either method
      plot_df$is_top <- (plot_df$rank_x <= highlight_top) |
        (plot_df$rank_y <= highlight_top)

      # Identify concordant top (top in both)
      plot_df$concordant_top <- (plot_df$rank_x <= highlight_top) &
        (plot_df$rank_y <= highlight_top)

      # Compute correlation for this pair
      rho <- cor(plot_df$x, plot_df$y, use = "complete.obs", method = cor_method)

      p <- ggplot(plot_df, aes(x = x, y = y)) +
        geom_point(aes(color = case_when(
          concordant_top ~ "Top in both",
          is_top ~ "Top in one",
          TRUE ~ "Other"
        )), alpha = 0.6, size = 1.5) +
        scale_color_manual(
          values = c("Top in both" = "#D62728",
                     "Top in one" = "#FF7F0E",
                     "Other" = "#1F77B4"),
          name = paste0("Top ", highlight_top)
        ) +
        geom_smooth(method = "lm", se = FALSE, color = "grey30",
                    linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = Inf, y = -Inf,
                 label = sprintf("ρ = %.3f", rho),
                 hjust = 1.1, vjust = -0.5, size = 4, fontface = "italic") +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5)) +
        labs(x = m1, y = m2,
             title = paste(m1, "vs", m2))

      return(p)
    })

    # Combine plots using patchwork if available, otherwise return list
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- patchwork::wrap_plots(plot_list, ncol = min(3, length(plot_list)))
      return(combined)
    } else {
      return(plot_list)
    }
  }

  # --- Generate requested plots ---
  if (plot_type == "matrix") {
    return(plot_matrix())
  } else if (plot_type == "scatter") {
    return(plot_scatter())
  } else {
    return(list(
      matrix = plot_matrix(),
      scatter = plot_scatter()
    ))
  }
}

#' Visualise method clustering based on correlations
#'
#' Creates a dendrogram showing hierarchical clustering of scoring methods
#' based on their pairwise correlations, with optional annotation of
#' modern vs. legacy classification.
#'
#' @param agreement Output from summarize_score_agreement()
#' @param highlight_groups Logical. Color branches by modern/legacy (default TRUE)
#'
#' @return A ggplot object (if ggdendro available) or base R plot
#' @export
plot_method_clustering <- function(agreement, highlight_groups = TRUE) {

  if (is.null(agreement$method_clusters)) {
    stop("No clustering information available. Need >= 3 methods.")
  }

  hc <- agreement$method_clusters
  meta <- agreement$metadata

  # Try ggdendro for prettier plots
  if (requireNamespace("ggdendro", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE)) {

    dend_data <- ggdendro::dendro_data(hc)

    # Create label data with method type
    labels_df <- dend_data$labels
    labels_df$method_type <- sapply(as.character(labels_df$label), function(m) {
      if (m %in% meta$modern_methods) return("Modern")
      if (m %in% meta$legacy_methods) return("Legacy")
      return("Other")
    })

    p <- ggplot2::ggplot() +
      geom_segment(data = dend_data$segments,
                             ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::geom_text(data = labels_df,
                         ggplot2::aes(x = x, y = y, label = label, color = method_type),
                         hjust = 1, angle = 90, size = 4, nudge_y = -0.02) +
      ggplot2::scale_color_manual(values = c("Modern" = "#2166AC",
                                             "Legacy" = "#B2182B",
                                             "Other" = "#666666"),
                                  name = "Method Type") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.2, 0.05))) +
      ggplot2::labs(title = "Scoring Method Clustering",
                    subtitle = "Based on rank correlation (1 - ρ)",
                    y = "Distance (1 - correlation)",
                    x = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    return(p)

  } else {
    # Fall back to base R
    plot(hc, main = "Scoring Method Clustering",
         sub = "Based on rank correlation (1 - ρ)",
         xlab = "", ylab = "Distance (1 - correlation)")

    invisible(hc)
  }
}


#' Create correlation heatmap with method type annotations
#'
#' Enhanced correlation matrix visualisation that highlights the
#' modern/legacy method structure.
#'
#' @param agreement Output from summarize_score_agreement()
#' @param order_by Character. How to order methods: "cluster" (default), "type", or "none"
#'
#' @return A ggplot object
#' @export
plot_correlation_heatmap <- function(agreement,
                                     order_by = c("cluster", "type", "none")) {

  order_by <- match.arg(order_by)

  cor_mat <- agreement$correlation_matrix
  meta <- agreement$metadata

  # Determine method order
  if (order_by == "cluster" && !is.null(agreement$method_clusters)) {
    method_order <- agreement$method_clusters$labels[agreement$method_clusters$order]
  } else if (order_by == "type") {
    method_order <- c(meta$legacy_methods, meta$modern_methods,
                      setdiff(meta$methods, c(meta$legacy_methods, meta$modern_methods)))
    method_order <- intersect(method_order, rownames(cor_mat))
  } else {
    method_order <- rownames(cor_mat)
  }

  # Reorder matrix
  cor_mat <- cor_mat[method_order, method_order]

  # Convert to long format
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Method1", "Method2", "Correlation")

  # Set factor levels for ordering
  cor_df$Method1 <- factor(cor_df$Method1, levels = method_order)
  cor_df$Method2 <- factor(cor_df$Method2, levels = rev(method_order))

  # Add method type for annotation
  cor_df$Type1 <- sapply(as.character(cor_df$Method1), function(m) {
    if (m %in% meta$modern_methods) return("Modern")
    if (m %in% meta$legacy_methods) return("Legacy")
    return("Other")
  })

  # Create text color based on correlation magnitude
  cor_df$text_color <- ifelse(abs(cor_df$Correlation) > 0.6, "white", "black")

  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Method1, y = Method2, fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Correlation)),
                       color = cor_df$text_color, size = 3.5) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-1, 1),
                                  name = paste0(tools::toTitleCase(meta$cor_method),
                                                "\nCorrelation")) +
    ggplot2::labs(title = "On-Target Score Correlations",
                  subtitle = paste0("Modern methods: ", paste(meta$modern_methods, collapse = ", "),
                                    "\nLegacy methods: ", paste(meta$legacy_methods, collapse = ", "))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
    ) +
    ggplot2::coord_fixed()

  return(p)
}
