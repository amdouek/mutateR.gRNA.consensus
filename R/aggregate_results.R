#' Aggregate batch checkpoint files into final result
#'
#' Combines all batch checkpoint files from a parallelized analysis run
#' into a single \code{mutateR_consensus_batch} object.
#'
#' @param output_dir Directory containing the checkpoints subdirectory.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return Path to the aggregated results file (invisibly).
#'
#' @details
#' This function is called automatically at the end of \code{\link{run_exome_analysis}},
#' but can also be called manually to re-aggregate results or after manually
#' adding/removing checkpoint files.
#'
#' @examples
#' \dontrun{
#' # Re-aggregate results after manual checkpoint manipulation
#' aggregate_batch_results("exome_analysis")
#'
#' # Load the aggregated results
#' results <- readRDS("exome_analysis/final_results.rds")
#' }
#'
#' @export
aggregate_batch_results <- function(output_dir, quiet = FALSE) {

  checkpoint_dir <- file.path(output_dir, "checkpoints")

  if (!dir.exists(checkpoint_dir)) {
    stop("Checkpoint directory not found: ", checkpoint_dir)
  }

  batch_files <- list.files(checkpoint_dir,
                            pattern = "^batch_.*\\.rds$",
                            full.names = TRUE)

  if (length(batch_files) == 0) {
    warning("No batch files found in ", checkpoint_dir)
    return(invisible(NULL))
  }

  # Sort batch files numerically
  batch_nums <- as.integer(gsub("batch_|\\.rds", "", basename(batch_files)))
  batch_files <- batch_files[order(batch_nums)]

  if (!quiet) {
    message("Aggregating ", length(batch_files), " batch files...")
  }

  # Load and combine all batches
  all_results <- list()

  for (bf in batch_files) {
    batch_data <- readRDS(bf)
    all_results <- c(all_results, batch_data)
  }

  if (!quiet) {
    message("Total genes aggregated: ", length(all_results))
  }

  # Load progress file for metadata
  progress_file <- file.path(output_dir, "progress.rds")
  progress <- if (file.exists(progress_file)) readRDS(progress_file) else NULL

  # Build batch result object
  final_result <- structure(
    all_results,
    class = c("mutateR_consensus_batch", "list"),
    metadata = list(
      gene_ids = names(all_results),
      n_successful = length(all_results),
      n_failed = if (!is.null(progress)) length(progress$failed_genes) else 0,
      failed_genes = if (!is.null(progress)) names(progress$failed_genes) else character(),
      species = if (!is.null(progress)) progress$parameters$species else NA_character_,
      nuclease = if (!is.null(progress)) progress$parameters$nuclease else NA_character_,
      methods = if (!is.null(progress)) progress$parameters$methods else NA_character_,
      modern_methods = if (!is.null(progress)) progress$parameters$modern_methods else NA_character_,
      legacy_methods = if (!is.null(progress)) progress$parameters$legacy_methods else NA_character_,
      cor_method = "spearman",
      analysis_type = "exome_scale",
      started_at = if (!is.null(progress)) progress$started_at else NA,
      completed_at = Sys.time()
    )
  )

  # Save final result
  final_path <- file.path(output_dir, "final_results.rds")
  saveRDS(final_result, final_path)

  if (!quiet) {
    message("Saved to: ", final_path)
  }

  return(invisible(final_path))
}


#' Extract aggregated scores from batch results
#'
#' Combines score data from all genes into a single data frame suitable
#' for SHAP analysis or other downstream processing.
#'
#' @param batch_result A \code{mutateR_consensus_batch} object, or path to RDS file.
#' @param include_sequences Logical. Include protospacer sequences (default TRUE).
#' @param include_coordinates Logical. Include genomic coordinates (default FALSE).
#'
#' @return A data frame with all gRNA scores across all genes.
#'
#' @examples
#' \dontrun{
#' results <- readRDS("exome_analysis/final_results.rds")
#' all_scores <- extract_all_scores(results)
#' dim(all_scores)  # Potentially millions of rows
#' }
#'
#' @export
extract_all_scores <- function(batch_result,
                               include_sequences = TRUE,
                               include_coordinates = FALSE) {

  # Load from file if path provided
  if (is.character(batch_result) && file.exists(batch_result)) {
    batch_result <- readRDS(batch_result)
  }

  if (!inherits(batch_result, "mutateR_consensus_batch")) {
    stop("Input must be a mutateR_consensus_batch object or path to one")
  }

  # Extract scores from each gene
  score_list <- lapply(names(batch_result), function(gene) {

    obj <- batch_result[[gene]]

    if (is.null(obj$scores) || nrow(obj$scores) == 0) {
      return(NULL)
    }

    df <- obj$scores

    # Add gene identifier
    df$gene <- gene
    df$gene_symbol <- obj$metadata$gene_symbol

    # Select columns
    keep_cols <- c("gene", "gene_symbol", "grna_id", "exon_rank")

    if (include_sequences) {
      keep_cols <- c(keep_cols, "protospacer_sequence", "sequence_context")
    }

    if (include_coordinates) {
      keep_cols <- c(keep_cols, "chr", "start", "end", "strand")
    }

    # Add method scores
    methods <- attr(df, "methods")
    if (is.null(methods)) {
      methods <- intersect(c("ruleset1", "deepspcas9", "ruleset3", "deephf",
                             "deepcpf1", "enpamgb"), names(df))
    }
    keep_cols <- c(keep_cols, methods)

    # Filter to available columns
    keep_cols <- intersect(keep_cols, names(df))

    df[, keep_cols, drop = FALSE]
  })

  # Remove NULLs and combine
  score_list <- Filter(Negate(is.null), score_list)

  if (length(score_list) == 0) {
    warning("No score data found in batch results")
    return(data.frame())
  }

  # Combine all data frames
  all_scores <- do.call(rbind, score_list)
  rownames(all_scores) <- NULL

  return(all_scores)
}


#' Compute summary statistics for batch results
#'
#' Generates per-gene summary statistics from exome-scale analysis results.
#'
#' @param batch_result A \code{mutateR_consensus_batch} object, or path to RDS file.
#'
#' @return A data frame with one row per gene and summary columns.
#'
#' @export
summarize_batch_results <- function(batch_result) {

  # Load from file if path provided
  if (is.character(batch_result) && file.exists(batch_result)) {
    batch_result <- readRDS(batch_result)
  }

  meta <- attr(batch_result, "metadata")
  methods <- if (!is.null(meta$methods)) meta$methods else
    c("ruleset1", "deepspcas9", "ruleset3", "deephf")

  # Summarize each gene
  summaries <- lapply(names(batch_result), function(gene) {

    obj <- batch_result[[gene]]

    base_info <- data.frame(
      gene = gene,
      gene_symbol = obj$metadata$gene_symbol %||% gene,
      transcript_id = obj$metadata$transcript_id %||% NA_character_,
      n_exons = if (!is.null(obj$exons)) length(obj$exons) else NA_integer_,
      n_sites = obj$metadata$n_sites %||% 0L,
      stringsAsFactors = FALSE
    )

    # Add agreement statistics if available
    if (!is.null(obj$agreement)) {
      cor_df <- obj$agreement$correlations

      if (!is.null(cor_df) && nrow(cor_df) > 0) {
        # Modern-only correlations
        modern_cors <- cor_df$correlation[cor_df$comparison_type == "modern-modern"]

        base_info$mean_cor_modern <- if (length(modern_cors) > 0) {
          mean(modern_cors, na.rm = TRUE)
        } else NA_real_

        base_info$min_cor_modern <- if (length(modern_cors) > 0) {
          min(modern_cors, na.rm = TRUE)
        } else NA_real_

        # Legacy-modern correlations
        legacy_cors <- cor_df$correlation[cor_df$comparison_type == "legacy-modern"]

        base_info$mean_cor_legacy <- if (length(legacy_cors) > 0) {
          mean(legacy_cors, na.rm = TRUE)
        } else NA_real_
      }
    }

    # Add per-method score statistics
    if (!is.null(obj$scores) && nrow(obj$scores) > 0) {
      for (m in methods) {
        if (m %in% names(obj$scores)) {
          scores <- obj$scores[[m]]
          base_info[[paste0("mean_", m)]] <- mean(scores, na.rm = TRUE)
          base_info[[paste0("median_", m)]] <- median(scores, na.rm = TRUE)
          base_info[[paste0("max_", m)]] <- max(scores, na.rm = TRUE)
        }
      }
    }

    return(base_info)
  })

  # Combine and return
  result <- do.call(rbind, summaries)
  rownames(result) <- NULL

  return(result)
}
