#' Print method for score agreement summary
#'
#' @param x Output from summarize_score_agreement()
#' @param ... Additional arguments (ignored)
#' @export
print.score_agreement <- function(x, ...) {

  meta <- x$metadata

  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║       mutateR On-Target Score Agreement Summary (Tiered)      ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

  # --- Dataset Overview ---
  cat("DATASET OVERVIEW\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("gRNAs analysed:    ", meta$n_grnas, "\n", sep = "")
  cat("Nuclease:          ", if (!is.null(meta$nuclease)) meta$nuclease else "Unknown", "\n", sep = "")
  cat("Correlation type:  ", meta$cor_method, "\n", sep = "")
  cat("\n")

  # --- Method Classification ---
  cat("METHOD CLASSIFICATION\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  if (length(meta$modern_methods) > 0) {
    cat("Modern methods:    ", paste(meta$modern_methods, collapse = ", "), "\n", sep = "")
  }
  if (length(meta$legacy_methods) > 0) {
    cat("Legacy methods:    ", paste(meta$legacy_methods, collapse = ", "), "\n", sep = "")
  }
  other_methods <- setdiff(meta$methods, c(meta$modern_methods, meta$legacy_methods))
  if (length(other_methods) > 0) {
    cat("Other methods:     ", paste(other_methods, collapse = ", "), "\n", sep = "")
  }
  cat("\n")

  # --- Tier Summary ---
  if (!is.null(x$tier_summary)) {
    cat("TIERED CORRELATION SUMMARY\n")
    cat("─────────────────────────────────────────────────────────────────\n")

    if (!is.null(x$tier_summary$modern)) {
      ts <- x$tier_summary$modern
      cat(sprintf("Modern methods (n=%d):   ρ = %.3f [%.3f - %.3f]\n",
                  length(ts$methods), ts$mean_correlation,
                  ts$min_correlation, ts$max_correlation))
    }

    if (!is.null(x$tier_summary$legacy_modern)) {
      ts <- x$tier_summary$legacy_modern
      cat(sprintf("Legacy vs Modern:        ρ = %.3f [%.3f - %.3f]\n",
                  ts$mean_correlation, ts$min_correlation, ts$max_correlation))
    }
    cat("\n")
  }

  # --- Pairwise Correlations ---
  cat("PAIRWISE CORRELATIONS\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  cor_display <- x$correlations[, c("method1", "method2", "correlation", "comparison_type")]
  cor_display$correlation <- sprintf("%.3f", cor_display$correlation)
  print(cor_display, row.names = FALSE)
  cat("\n")

  # --- Top-N Overlap ---
  cat("TOP-N RANK OVERLAP\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  for (n in unique(x$rank_overlaps$top_n)) {
    sub <- x$rank_overlaps[x$rank_overlaps$top_n == n, ]

    # Overall
    mean_enrich_all <- mean(sub$enrichment, na.rm = TRUE)

    # Modern-only
    modern_sub <- sub[sub$comparison_type == "modern-modern", ]
    mean_enrich_modern <- if (nrow(modern_sub) > 0) {
      mean(modern_sub$enrichment, na.rm = TRUE)
    } else NA

    cat(sprintf("  Top %d:\n", n))
    cat(sprintf("    All pairs:     %.1fx enrichment (vs. random)\n", mean_enrich_all))
    if (!is.na(mean_enrich_modern)) {
      cat(sprintf("    Modern only:   %.1fx enrichment\n", mean_enrich_modern))
    }
  }
  cat("\n")

  # --- Concordant High Performers ---
  cat("CONCORDANT HIGH PERFORMERS (Top 5)\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("(gRNAs ranking well across methods - lower max_rank = more consistent)\n\n")

  if (!is.null(x$concordant_top) && nrow(x$concordant_top) >= 1) {
    top5 <- head(x$concordant_top, 5)
    # Select display columns
    display_cols <- c("grna_id", "protospacer_sequence")
    rank_cols <- grep("^rank_", names(top5), value = TRUE)
    max_col <- grep("max_rank", names(top5), value = TRUE)[1]
    display_cols <- c(display_cols, rank_cols, max_col)
    display_cols <- intersect(display_cols, names(top5))
    print(top5[, display_cols], row.names = FALSE)
  } else {
    cat("(No data available)\n")
  }
  cat("\n")

  # --- Modern Method Discordance ---
  if (!is.null(x$discordant_examples) && nrow(x$discordant_examples) >= 1) {
    cat("MOST DISCORDANT - MODERN METHODS (Top 5)\n")
    cat("─────────────────────────────────────────────────────────────────\n")
    cat("(gRNAs with largest rank disagreement among modern methods)\n\n")

    disc5 <- head(x$discordant_examples, 5)
    display_cols <- c("grna_id", "protospacer_sequence")
    rank_cols <- grep("^rank_", names(disc5), value = TRUE)
    display_cols <- c(display_cols, rank_cols, "rank_variance_modern")
    display_cols <- intersect(display_cols, names(disc5))
    print(disc5[, display_cols], row.names = FALSE)
    cat("\n")
  }

  # --- Legacy Divergence ---
  if (!is.null(x$legacy_divergent) && nrow(x$legacy_divergent) >= 1) {
    cat("LEGACY METHOD DIVERGENCE (Top 5)\
")
    cat("─────────────────────────────────────────────────────────────────\n")
    cat("(gRNAs where legacy method most disagrees with modern consensus)\n")
    cat("(Negative divergence = legacy ranks it HIGHER than modern methods)\n\n")

    leg5 <- head(x$legacy_divergent, 5)
    print(leg5, row.names = FALSE)
    cat("\n")
  }

  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Access tiered metrics via: extract_tiered_metrics(x, tier = 'modern')\n")
  cat("Available tiers: 'modern', 'robust', 'all', 'legacy'\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  invisible(x)
}


#' Summary method for score agreement
#' @export
summary.score_agreement <- function(object, ...) {

  meta <- object$metadata

  cat("\nScore Agreement Summary\n")
  cat("=======================\n\n")

  cat("Methods: ", paste(meta$methods, collapse = ", "), "\n")
  cat("gRNAs:   ", meta$n_grnas, "\n\n")

  # Correlation summary by type
  cat("Correlations by comparison type:\n")
  cor_by_type <- aggregate(correlation ~ comparison_type,
                           data = object$correlations,
                           FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  print(cor_by_type)

  invisible(object)
}
