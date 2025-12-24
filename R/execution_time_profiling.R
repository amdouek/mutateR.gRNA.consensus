#' Profile execution time for consensus analysis steps
#'
#' Runs a single gene analysis with timing for each step to identify bottlenecks.
#'
#' @param gene_id Gene symbol or Ensembl ID
#' @param species Species code
#' @param genome BSgenome object
#' @param nuclease Cas effector
#' @param methods Scoring methods
#'
#' @return List containing timing breakdown and results
#' @export
profile_analysis_timing <- function(gene_id,
                                    species,
                                    genome,
                                    nuclease = "Cas9",
                                    methods = c("ruleset1", "deepspcas9", "ruleset3", "deephf")) {

  timings <- list()

  # --- Step 1: Gene Info ---
  timings$gene_info <- system.time({
    tx_info <- get_gene_info(gene_id, species, id_type = "symbol")
    tx_id <- tx_info$canonical$ensembl_transcript_id[1]
  })["elapsed"]

  # --- Step 2: Exon Structures ---
  timings$exon_structures <- system.time({
    exons_gr <- get_exon_structures(tx_id, species, output = "GRanges")
  })["elapsed"]

  # --- Step 3: Find Sites ---
  timings$find_sites <- system.time({
    sites_gr <- find_cas9_sites(exons_gr, genome, pam = "NGG", protospacer_length = 20)
  })["elapsed"]

  n_sites <- length(sites_gr)

  # --- Step 4: Individual Scoring Methods ---
  timings$scoring <- list()

  for (method in methods) {
    timings$scoring[[method]] <- system.time({
      scored <- score_grnas(sites_gr, method = method)
    })["elapsed"]
  }

  # --- Step 5: Agreement Analysis ---
  # First get multi-scores
  timings$multi_score_combine <- system.time({
    multi_scores <- score_grnas_multi(sites_gr, methods = methods, quiet = TRUE)
  })["elapsed"]

  timings$agreement_analysis <- system.time({
    agreement <- summarize_score_agreement(multi_scores)
  })["elapsed"]

  # --- Step 6: Plotting ---
  timings$plotting <- system.time({
    p1 <- plot_score_correlations(multi_scores, plot_type = "matrix")
    p2 <- plot_correlation_heatmap(agreement)
  })["elapsed"]

  # --- Summary ---
  total_scoring <- sum(unlist(timings$scoring))
  total_time <- sum(unlist(timings[c("gene_info", "exon_structures", "find_sites",
                                     "agreement_analysis", "plotting")])) + total_scoring

  summary_df <- data.frame(
    step = c("Gene Info (API)", "Exon Structures (API)", "Find Sites",
             paste0("Score: ", methods), "Agreement Analysis", "Plotting", "TOTAL"),
    seconds = c(timings$gene_info, timings$exon_structures, timings$find_sites,
                unlist(timings$scoring), timings$agreement_analysis, timings$plotting, total_time),
    stringsAsFactors = FALSE
  )

  summary_df$percent <- round(100 * summary_df$seconds / total_time, 1)

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Timing Profile for ", gene_id, " (", n_sites, " gRNAs)\n", sep = "")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  print(summary_df, row.names = FALSE)
  cat("\n")

  # Identify bottleneck
  bottleneck_idx <- which.max(summary_df$seconds[-nrow(summary_df)])
  cat("Primary bottleneck: ", summary_df$step[bottleneck_idx],
      " (", summary_df$percent[bottleneck_idx], "%)\n\n", sep = "")

  return(list(
    timings = timings,
    summary = summary_df,
    n_sites = n_sites,
    bottleneck = summary_df$step[bottleneck_idx]
  ))
}
