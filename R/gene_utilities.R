#' Get all protein-coding genes for a species
#'
#' Retrieves gene symbols for all protein-coding genes via BioMart.
#' Results are cached locally to avoid repeated API calls (reducing runtime).
#'
#' @param species Character. Ensembl species code (default "hsapiens").
#' @param cache Logical. Cache results locally (default TRUE).
#' @param cache_dir Character. Cache directory (default "~/.mutateR_cache").
#' @param max_age_days Numeric. Maximum age of cache before refresh (default 30).
#'
#' @return Character vector of gene symbols.
#'
#' @examples
#' \dontrun{
#' # Get all human protein-coding genes
#' human_genes <- get_protein_coding_genes("hsapiens")
#' length(human_genes)  # ~20,000
#'
#' # Get mouse genes
#' mouse_genes <- get_protein_coding_genes("mmusculus")
#' }
#'
#' @export
get_protein_coding_genes <- function(species = "hsapiens",
                                     cache = TRUE,
                                     cache_dir = "~/.mutateR_cache",
                                     max_age_days = 30) {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' required. Install with: BiocManager::install('biomaRt')")
  }

  # Expand path
  cache_dir <- path.expand(cache_dir)

  # Check cache
  if (cache) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, paste0(species, "_protein_coding_genes.rds"))

    if (file.exists(cache_file)) {
      file_age <- difftime(Sys.time(), file.info(cache_file)$mtime, units = "days")

      if (file_age < max_age_days) {
        message("Loading cached gene list (", round(file_age, 1), " days old)")
        return(readRDS(cache_file))
      } else {
        message("Cache expired (", round(file_age, 1), " days old). Refreshing...")
      }
    }
  }

  message("Fetching protein-coding genes from Ensembl BioMart...")

  # Connect to Ensembl
  ensembl <- tryCatch({
    biomaRt::useEnsembl(
      biomart = "genes",
      dataset = paste0(species, "_gene_ensembl")
    )
  }, error = function(e) {
    stop("Failed to connect to Ensembl BioMart: ", e$message,
         "\nCheck species code and internet connection.")
  })

  # Query protein-coding genes
  genes <- biomaRt::getBM(
    attributes = c("external_gene_name", "ensembl_gene_id", "gene_biotype"),
    filters = "biotype",
    values = "protein_coding",
    mart = ensembl
  )

  # Extract unique gene symbols (excluding empty/NA)
  gene_symbols <- unique(genes$external_gene_name)
  gene_symbols <- gene_symbols[gene_symbols != "" & !is.na(gene_symbols)]

  # Sort alphabetically
  gene_symbols <- sort(gene_symbols)

  message("Found ", length(gene_symbols), " protein-coding genes")

  # Cache results
  if (cache) {
    saveRDS(gene_symbols, cache_file)
    message("Cached to: ", cache_file)
  }

  return(gene_symbols)
}


#' Monitor progress of running exome analysis
#'
#' Displays current status of an in-progress or completed exome analysis,
#' including completion percentage, rate, and ETA.
#'
#' @param output_dir Analysis output directory.
#'
#' @return The progress object (invisibly).
#'
#' @examples
#' \dontrun{
#' # From a separate R session while analysis is running:
#' monitor_exome_analysis("exome_analysis")
#' }
#'
#' @export
monitor_exome_analysis <- function(output_dir) {

  progress_file <- file.path(output_dir, "progress.rds")

  if (!file.exists(progress_file)) {
    message("No progress file found at: ", progress_file)
    message("Analysis may not have started yet.")
    return(invisible(NULL))
  }

  progress <- readRDS(progress_file)

  elapsed <- as.numeric(difftime(Sys.time(), progress$started_at, units = "hours"))
  completed <- length(progress$completed_genes)
  failed <- length(progress$failed_genes)
  total <- progress$total_genes
  remaining <- total - completed - failed

  rate <- if (elapsed > 0) completed / elapsed else 0
  eta_hours <- if (rate > 0) remaining / rate else NA

  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║           Exome Analysis Progress Monitor                     ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

  cat("TIMING\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("Started:          ", format(progress$started_at, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Elapsed:          ", sprintf("%.1f hours", elapsed), "\n")
  cat("Last checkpoint:  ", format(progress$last_checkpoint, "%Y-%m-%d %H:%M:%S"), "\n")

  # Check if analysis appears stalled
  since_checkpoint <- as.numeric(difftime(Sys.time(), progress$last_checkpoint, units = "mins"))
  if (since_checkpoint > 30 && remaining > 0) {
    cat("  ⚠ No checkpoint in ", round(since_checkpoint), " minutes - analysis may be stalled\n")
  }
  cat("\n")

  cat("PROGRESS\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("Total genes:      ", total, "\n")
  cat("Completed:        ", completed, sprintf(" (%.1f%%)", 100 * completed / total), "\n")
  cat("Failed:           ", failed, "\n")
  cat("Remaining:        ", remaining, "\n\n")

  cat("PERFORMANCE\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("Rate:             ", sprintf("%.1f genes/hour", rate), "\n")

  if (!is.na(eta_hours) && remaining > 0) {
    eta_time <- Sys.time() + eta_hours * 3600
    cat("ETA:              ", format(eta_time, "%Y-%m-%d %H:%M"),
        sprintf(" (%.1f hours)", eta_hours), "\n")
  } else if (remaining == 0) {
    cat("Status:           Complete!\n")
  }
  cat("\n")

  # Progress bar
  pct <- completed / total
  bar_width <- 50
  filled <- round(pct * bar_width)

  bar <- paste0(
    "[",
    paste(rep("█", filled), collapse = ""),
    paste(rep("░", bar_width - filled), collapse = ""),
    "]"
  )

  cat("PROGRESS BAR\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(bar, sprintf(" %.1f%%\n", 100 * pct))

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  return(invisible(progress))
}


#' Retry failed genes from a previous run
#'
#' Attempts to reprocess genes that failed during a previous analysis run.
#' Useful for recovering from transient errors (e.g., API timeouts).
#' You should check the failed gene list to see if these failure cases are due
#' to the gene or the pipeline run.
#'
#' @param output_dir Analysis output directory.
#' @param genome BSgenome object.
#' @param n_workers Number of workers (default 4, lower than initial run).
#' @param clear_failures Logical. Clear previous failure records before retry (default TRUE).
#'
#' @return Number of genes successfully recovered.
#'
#' @examples
#' \dontrun{
#' # After initial run completes with some failures
#' recovered <- retry_failed_genes("exome_analysis", BSgenome.Hsapiens.UCSC.hg38)
#' message(recovered, " genes recovered on retry")
#' }
#'
#' @export
retry_failed_genes <- function(output_dir,
                               genome,
                               n_workers = 4,
                               clear_failures = TRUE) {

  progress_file <- file.path(output_dir, "progress.rds")

  if (!file.exists(progress_file)) {
    stop("No progress file found in ", output_dir)
  }

  progress <- readRDS(progress_file)
  failed_genes <- names(progress$failed_genes)

  if (length(failed_genes) == 0) {
    message("No failed genes to retry!")
    return(invisible(0))
  }

  message("Found ", length(failed_genes), " previously failed genes")

  # Show failure reasons
  message("\nFailure summary:")
  failure_reasons <- unlist(progress$failed_genes)
  reason_table <- table(failure_reasons)
  for (reason in names(head(sort(reason_table, decreasing = TRUE), 5))) {
    message("  ", reason_table[reason], "x: ", substr(reason, 1, 60))
  }
  message("")

  # Clear failed genes from progress if requested
  if (clear_failures) {
    progress$failed_genes <- list()
    saveRDS(progress, progress_file)
  }

  # Re-run with failed genes
  run_exome_analysis(
    gene_ids = failed_genes,
    species = progress$parameters$species,
    genome = genome,
    output_dir = output_dir,
    n_workers = n_workers,
    batch_size = min(50, length(failed_genes)),
    resume = TRUE,
    nuclease = progress$parameters$nuclease,
    methods = progress$parameters$methods,
    skip_plots = TRUE
  )

  # Check how many were recovered
  progress_new <- readRDS(progress_file)
  still_failed <- length(progress_new$failed_genes)
  recovered <- length(failed_genes) - still_failed

  message("\n")
  message("═══════════════════════════════════════════════════════════════")
  message("Retry Results")
  message("═══════════════════════════════════════════════════════════════")
  message("Attempted:  ", length(failed_genes))
  message("Recovered:  ", recovered)
  message("Still failed: ", still_failed)

  if (still_failed > 0) {
    message("\nPersistently failing genes saved to: ",
            file.path(output_dir, "failed_genes.rds"))
  }

  return(invisible(recovered))
}


#' Get analysis completion status
#'
#' Quick check of whether an analysis is complete, in progress, or not started.
#'
#' @param output_dir Analysis output directory.
#'
#' @return Character: "complete", "in_progress", "not_started", or "error".
#'
#' @export
get_analysis_status <- function(output_dir) {

  progress_file <- file.path(output_dir, "progress.rds")
  final_file <- file.path(output_dir, "final_results.rds")

  if (!file.exists(progress_file)) {
    return("not_started")
  }

  if (file.exists(final_file)) {
    return("complete")
  }

  # Check if in progress
  progress <- tryCatch({
    readRDS(progress_file)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(progress)) {
    return("error")
  }

  completed <- length(progress$completed_genes)
  failed <- length(progress$failed_genes)
  total <- progress$total_genes

  if (completed + failed >= total) {
    return("complete")
  }

  return("in_progress")
}


#' List available cached gene lists
#'
#' Shows what species have cached protein-coding gene lists available.
#'
#' @param cache_dir Character. Cache directory (default "~/.mutateR_cache").
#'
#' @return Data frame with species and cache information.
#'
#' @export
list_cached_genes <- function(cache_dir = "~/.mutateR_cache") {

  cache_dir <- path.expand(cache_dir)

  if (!dir.exists(cache_dir)) {
    message("Cache directory does not exist: ", cache_dir)
    return(data.frame())
  }

  cache_files <- list.files(cache_dir,
                            pattern = "_protein_coding_genes\\.rds$",
                            full.names = TRUE)

  if (length(cache_files) == 0) {
    message("No cached gene lists found")
    return(data.frame())
  }

  # Build summary
  summaries <- lapply(cache_files, function(f) {
    species <- gsub("_protein_coding_genes\\.rds$", "", basename(f))
    info <- file.info(f)
    genes <- readRDS(f)

    data.frame(
      species = species,
      n_genes = length(genes),
      cache_date = as.Date(info$mtime),
      age_days = as.numeric(difftime(Sys.time(), info$mtime, units = "days")),
      size_mb = round(info$size / 1024^2, 2),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, summaries)
  rownames(result) <- NULL

  return(result)
}


#' Null-coalescing operator
#'
#' Returns the left-hand side if not NULL, otherwise the right-hand side.
#'
#' @param x Left-hand side value.
#' @param y Right-hand side (default) value.
#'
#' @return x if not NULL, otherwise y.
#'
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
