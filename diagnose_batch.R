#' Diagnose batch analysis completeness
#'
#' @param batch_result A mutateR_consensus_batch object
#' @return Data frame summarising each gene's status
#' @export
diagnose_batch_result <- function(batch_result) {

  meta <- attr(batch_result, "metadata")

  if (is.null(meta)) {
    stop("Input does not appear to be a mutateR_consensus_batch object")
  }

  # Check explicitly failed genes
  failed <- meta$failed_genes

  # Check genes in results
  successful <- names(batch_result)

  # For each successful gene, check data completeness
  status_list <- lapply(successful, function(gene) {
    obj <- batch_result[[gene]]

    data.frame(
      gene = gene,
      status = "success",
      n_exons = if (!is.null(obj$exons)) length(obj$exons) else NA_integer_,
      n_sites = obj$metadata$n_sites,
      has_scores = !is.null(obj$scores) && nrow(obj$scores) > 0,
      n_scored = if (!is.null(obj$scores)) nrow(obj$scores) else 0L,
      has_agreement = !is.null(obj$agreement),
      stringsAsFactors = FALSE
    )
  })

  status_df <- do.call(rbind, status_list)

  # Add failed genes
  if (length(failed) > 0) {
    failed_df <- data.frame(
      gene = failed,
      status = "failed",
      n_exons = NA_integer_,
      n_sites = NA_integer_,
      has_scores = FALSE,
      n_scored = 0L,
      has_agreement = FALSE,
      stringsAsFactors = FALSE
    )
    status_df <- rbind(status_df, failed_df)
  }

  # Check for genes with 0 scores (may have succeeded but with no data)
  status_df$usable_for_shap <- status_df$n_scored > 0

  # Summary
  cat("Batch Diagnosis\n")
  cat("===============\n")
  cat("Requested:", meta$n_requested, "genes\n")
  cat("Successful:", meta$n_successful, "genes\n")
  cat("Failed:", meta$n_failed, "genes\n")
  cat("Usable for SHAP:", sum(status_df$usable_for_shap), "genes\n")

  if (meta$n_failed > 0) {
    cat("\nFailed genes:", paste(failed, collapse = ", "), "\n")
  }

  # Check for zero-site genes
  zero_site <- status_df$gene[status_df$n_scored == 0 & status_df$status == "success"]
  if (length(zero_site) > 0) {
    cat("\nGenes with 0 gRNA sites:", paste(zero_site, collapse = ", "), "\n")
  }

  cat("\n")

  return(status_df)
}
