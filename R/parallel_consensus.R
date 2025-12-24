#' Run exome-scale parallel consensus analysis
#'
#' Processes thousands of genes in parallel with checkpointing for fault tolerance.
#' Designed for overnight/weekend runs on local hardware. Results are saved in
#' batches to enable resumption from interruptions.
#'
#' Note: Currently only focusing on the 'exome' per the current mutateR architecture
#' which currently only handles exonic sequences. A future version may allow for intronic targeting.
#'
#' @param gene_ids Character vector. Gene symbols to analyse.
#' @param species Character. Ensembl species code (default "hsapiens").
#' @param genome BSgenome object for the target species.
#' @param output_dir Character. Directory for checkpoints and results (default "exome_analysis").
#' @param n_workers Integer. Number of parallel workers (default 8).
#' @param batch_size Integer. Genes per checkpoint batch (default 100).
#' @param resume Logical. Resume from existing progress if available (default TRUE).
#' @param nuclease Character. One of "Cas9" (default), "Cas12a", or "enCas12a".
#' @param methods Character vector. Scoring methods to use. If NULL (default), auto-selects
#'        based on nuclease type.
#' @param skip_plots Logical. Skip plot generation for speed (default TRUE).
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return Path to final aggregated results file (invisibly).
#'
#' @details
#' The function processes genes in parallel batches, saving checkpoints after each
#' batch completes. If interrupted, subsequent runs with \code{resume = TRUE} will
#' continue from the last checkpoint without reprocessing completed genes.
#'
#' Recommended worker counts based on available RAM:
#' \itemize{
#'   \item 16GB RAM: 4-6 workers
#'   \item 32GB RAM: 8-10 workers
#'   \item 64GB RAM: 12-16 workers
#' }
#'
#' @section Output Structure:
#' The function creates the following directory structure:
#' \preformatted{
#' output_dir/
#' ├── checkpoints/
#' │   ├── batch_0001.rds
#' │   ├── batch_0002.rds
#' │   └── ...
#' ├── progress.rds
#' ├── failed_genes.rds
#' └── final_results.rds
#' }
#'
#' @examples
#' \dontrun{
#' library(mutateR)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Get protein-coding genes
#' all_genes <- get_protein_coding_genes("hsapiens")
#'
#' # Run full exome analysis (overnight)
#' run_exome_analysis(
#'   gene_ids = all_genes,
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   output_dir = "exome_analysis",
#'   n_workers = 8
#' )
#'
#' # Monitor progress from another session
#' monitor_exome_analysis("exome_analysis")
#'
#' # After completion, load results
#' results <- readRDS("exome_analysis/final_results.rds")
#' }
#'
#' @seealso
#' \code{\link{monitor_exome_analysis}} for progress monitoring,
#' \code{\link{retry_failed_genes}} for recovering failed genes,
#' \code{\link{aggregate_batch_results}} for manual result aggregation.
#'
#' @export
run_exome_analysis <- function(gene_ids,
                               species = "hsapiens",
                               genome,
                               output_dir = "exome_analysis",
                               n_workers = 8,
                               batch_size = 100,
                               resume = TRUE,
                               nuclease = c("Cas9", "Cas12a", "enCas12a"),
                               methods = NULL,
                               skip_plots = TRUE,
                               quiet = FALSE) {

  # ═══════════════════════════════════════════════════════════════
  # Validation and Setup
  # ═══════════════════════════════════════════════════════════════

  nuclease <- match.arg(nuclease)

  # Check required packages
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' required. Install with: install.packages('future')")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' required. Install with: install.packages('future.apply')")
  }

  # Validate genome
  if (!inherits(genome, "BSgenome")) {
    stop("Please supply a valid BSgenome object.")
  }

  # Validate gene_ids
  if (length(gene_ids) == 0) {
    stop("gene_ids cannot be empty.")
  }

  # Remove any duplicates
  gene_ids <- unique(gene_ids)

  # Set default methods based on nuclease
  if (is.null(methods)) {
    methods <- switch(nuclease,
                      "Cas9" = c("ruleset1", "deepspcas9", "ruleset3", "deephf"),
                      "Cas12a" = c("deepcpf1", "enpamgb"),
                      "enCas12a" = c("deepcpf1", "enpamgb"))
  }

  # Set method classification
  if (nuclease == "Cas9") {
    modern_methods <- intersect(c("deepspcas9", "deephf", "ruleset3"), methods)
    legacy_methods <- intersect(c("ruleset1"), methods)
  } else {
    modern_methods <- methods
    legacy_methods <- character(0)
  }

  # Create output directories
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  dir.create(checkpoint_dir, showWarnings = FALSE)

  progress_file <- file.path(output_dir, "progress.rds")
  failed_file <- file.path(output_dir, "failed_genes.rds")

  # ═══════════════════════════════════════════════════════════════
  # Progress Management
  # ═══════════════════════════════════════════════════════════════

  if (resume && file.exists(progress_file)) {
    progress <- readRDS(progress_file)

    if (!quiet) {
      message("Resuming from checkpoint...")
      message("  Previously completed: ", length(progress$completed_genes), " genes")
      message("  Previously failed: ", length(progress$failed_genes), " genes")
    }

    # Validate that parameters match
    if (!is.null(progress$parameters)) {
      if (progress$parameters$species != species) {
        warning("Species mismatch with previous run. Previous: ",
                progress$parameters$species, ", Current: ", species)
      }
      if (progress$parameters$nuclease != nuclease) {
        warning("Nuclease mismatch with previous run. Previous: ",
                progress$parameters$nuclease, ", Current: ", nuclease)
      }
    }

  } else {
    progress <- list(
      started_at = Sys.time(),
      total_genes = length(gene_ids),
      completed_genes = character(),
      completed_batches = integer(),
      failed_genes = list(),
      parameters = list(
        species = species,
        nuclease = nuclease,
        methods = methods,
        modern_methods = modern_methods,
        legacy_methods = legacy_methods,
        batch_size = batch_size,
        n_workers = n_workers
      ),
      last_checkpoint = Sys.time()
    )
  }

  # Identify remaining genes (deduplication)
  remaining_genes <- setdiff(gene_ids, progress$completed_genes)
  remaining_genes <- setdiff(remaining_genes, names(progress$failed_genes))

  if (length(remaining_genes) == 0) {
    if (!quiet) message("All genes already processed!")
    final_path <- aggregate_batch_results(output_dir, quiet = quiet)
    return(invisible(final_path))
  }

  # ═══════════════════════════════════════════════════════════════
  # Print Analysis Plan
  # ═══════════════════════════════════════════════════════════════

  if (!quiet) {
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════════╗\n")
    cat("║           mutateR Exome-Scale Analysis                        ║\n")
    cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

    cat("ANALYSIS CONFIGURATION\n")
    cat("─────────────────────────────────────────────────────────────────\n")
    cat("Species:                ", species, "\n")
    cat("Nuclease:               ", nuclease, "\n")
    cat("Methods:                ", paste(methods, collapse = ", "), "\n")
    cat("Output directory:       ", output_dir, "\n")
    cat("Workers:                ", n_workers, "\n")
    cat("Batch size:             ", batch_size, "\n")
    cat("Skip plots:             ", skip_plots, "\n\n")

    cat("PROGRESS\n")
    cat("─────────────────────────────────────────────────────────────────\n")
    cat("Total genes requested:  ", length(gene_ids), "\n")
    cat("Already completed:      ", length(progress$completed_genes), "\n")
    cat("Previously failed:      ", length(progress$failed_genes), "\n")
    cat("Remaining to process:   ", length(remaining_genes), "\n\n")

    # Time estimate
    estimated_seconds <- length(remaining_genes) * 37.5 / n_workers
    estimated_hours <- estimated_seconds / 3600

    cat("ESTIMATED TIME\n")
    cat("─────────────────────────────────────────────────────────────────\n")
    cat("Estimated time remaining: ", round(estimated_hours, 1), " hours\n")

    if (estimated_hours > 1) {
      eta <- Sys.time() + estimated_seconds
      cat("Estimated completion:     ", format(eta, "%Y-%m-%d %H:%M"), "\n")
    }
    cat("\n")
  }

  # ═══════════════════════════════════════════════════════════════
  # Batch Processing
  # ═══════════════════════════════════════════════════════════════

  # Split remaining genes into batches
  n_remaining <- length(remaining_genes)
  batch_indices <- split(seq_len(n_remaining),
                         ceiling(seq_len(n_remaining) / batch_size))
  n_batches <- length(batch_indices)

  # Determine batch numbering (continue from previous)
  start_batch_num <- max(c(0, progress$completed_batches)) + 1

  if (!quiet) {
    cat("BATCH PROCESSING\n")
    cat("─────────────────────────────────────────────────────────────────\n")
    cat("Processing ", n_batches, " batches (",
        start_batch_num, " to ", start_batch_num + n_batches - 1, ")\n\n", sep = "")
  }

  # ═══════════════════════════════════════════════════════════════
  # Python Environment Setup (CRITICAL FOR PARALLEL WORKERS)
  # ═══════════════════════════════════════════════════════════════
  #
  # Workers are fresh R sessions. If reticulate initializes before
  # we can activate the correct environment, it will use the wrong
  # Python and cannot be switched. We solve this by:
  # 1. Setting RETICULATE_PYTHON env var (inherited by workers)
  # 2. Setting RETICULATE_PYTHON_ENV env var as backup
  # 3. Workers will use these on their first reticulate call

  if (!quiet) {
    message("Configuring Python environment for parallel workers...")
  }

  # Activate in parent session first
  if (!mutateR::check_mutater_env()) {
    env_activated <- mutateR::activate_mutater_env()
    if (!env_activated) {
      stop("Could not activate mutateR Python environment.\n",
           "Please run mutateR::install_mutater_env() first.")
    }
  }

  # Get the Python executable path from the activated environment
  py_config <- reticulate::py_config()
  python_path <- py_config$python

  if (is.null(python_path) || !file.exists(python_path)) {
    stop("Could not determine Python path from activated environment.")
  }

  # Verify this is the r-mutater environment (not base miniconda)
  if (!grepl("r-mutater", python_path, fixed = TRUE)) {
    stop("Python path does not appear to be the r-mutater environment: ", python_path)
  }

  # Set environment variables - these propagate to worker processes
  Sys.setenv(RETICULATE_PYTHON = python_path)
  Sys.setenv(RETICULATE_PYTHON_ENV = "r-mutater")

  if (!quiet) {
    message("  RETICULATE_PYTHON set to: ", python_path)
  }

  # Verify TensorFlow is available
  if (!reticulate::py_module_available("tensorflow")) {
    stop("TensorFlow not found. Please run mutateR::install_mutater_env(fresh = TRUE)")
  }

  if (!quiet) {
    message("  TensorFlow: OK")
    message("")
  }

  # Setup parallel backend
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = n_workers)

  # Process each batch
  for (b in seq_along(batch_indices)) {
    batch_num <- start_batch_num + b - 1
    batch_genes <- remaining_genes[batch_indices[[b]]]

    batch_start_time <- Sys.time()

    if (!quiet) {
      cat("─────────────────────────────────────────────────────────────────\n")
      cat("Batch ", batch_num, "/", start_batch_num + n_batches - 1,
          " | ", length(batch_genes), " genes | Started: ",
          format(batch_start_time, "%H:%M:%S"), "\n", sep = "")
      cat("─────────────────────────────────────────────────────────────────\n")
    }

    # Process batch in parallel
    batch_results <- process_batch_parallel(
      gene_ids = batch_genes,
      species = species,
      genome = genome,
      nuclease = nuclease,
      methods = methods,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      skip_plots = skip_plots
    )

    # Separate successes and failures
    is_success <- sapply(batch_results, function(x) !isTRUE(x$error))

    successful_results <- batch_results[is_success]
    failed_in_batch <- batch_results[!is_success]

    # Update progress
    progress$completed_genes <- c(progress$completed_genes, names(successful_results))
    progress$completed_batches <- c(progress$completed_batches, batch_num)

    for (gene in names(failed_in_batch)) {
      progress$failed_genes[[gene]] <- failed_in_batch[[gene]]$message
    }

    progress$last_checkpoint <- Sys.time()

    # Save checkpoint
    if (length(successful_results) > 0) {
      batch_file <- file.path(checkpoint_dir, sprintf("batch_%04d.rds", batch_num))
      saveRDS(successful_results, batch_file)
    }

    saveRDS(progress, progress_file)

    # Batch timing report
    batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))
    genes_per_min <- length(batch_genes) / (batch_elapsed / 60)

    if (!quiet) {
      cat("  Completed:  ", length(successful_results), " genes\n", sep = "")
      cat("  Failed:     ", length(failed_in_batch), " genes\n", sep = "")
      cat("  Batch time: ", round(batch_elapsed / 60, 1), " minutes\n", sep = "")
      cat("  Rate:       ", round(genes_per_min, 1), " genes/minute\n", sep = "")

      # Overall progress
      total_done <- length(progress$completed_genes)
      total_pct <- round(100 * total_done / length(gene_ids), 1)
      cat("  Progress:   ", total_done, "/", length(gene_ids),
          " (", total_pct, "%)\n", sep = "")

      # ETA
      if (total_done > 0) {
        elapsed_total <- as.numeric(difftime(Sys.time(), progress$started_at, units = "secs"))
        rate <- total_done / elapsed_total
        remaining_count <- length(gene_ids) - total_done - length(progress$failed_genes)
        eta_seconds <- remaining_count / rate
        eta_time <- Sys.time() + eta_seconds
        cat("  ETA:        ", format(eta_time, "%Y-%m-%d %H:%M"), "\n", sep = "")
      }
      cat("\n")
    }

    # Force garbage collection between batches
    gc(verbose = FALSE)
  }

  # ═══════════════════════════════════════════════════════════════
  # Finalisation
  # ═══════════════════════════════════════════════════════════════

  total_elapsed <- as.numeric(difftime(Sys.time(), progress$started_at, units = "hours"))

  if (!quiet) {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("                    Analysis Complete!                         \n")
    cat("═══════════════════════════════════════════════════════════════\n\n")

    cat("Total genes processed: ", length(progress$completed_genes), "\n")
    cat("Total genes failed:    ", length(progress$failed_genes), "\n")
    cat("Total time:            ", round(total_elapsed, 2), " hours\n")
    cat("Average rate:          ", round(length(progress$completed_genes) / total_elapsed, 1),
        " genes/hour\n\n")
  }

  # Save failed genes summary
  if (length(progress$failed_genes) > 0) {
    saveRDS(progress$failed_genes, failed_file)
    if (!quiet) {
      cat("Failed genes saved to: ", failed_file, "\n\n")
    }
  }

  # Aggregate results
  if (!quiet) cat("Aggregating batch results...\n")
  final_path <- aggregate_batch_results(output_dir, quiet = quiet)

  if (!quiet) {
    cat("\nFinal results saved to: ", final_path, "\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")
  }

  return(invisible(final_path))
}


#' Process a batch of genes in parallel
#'
#' Internal function that processes a batch of genes using parallel workers.
#' Each gene is analysed independently, with errors captured rather than
#' halting the batch.
#'
#' @param gene_ids Character vector of gene symbols.
#' @param species Ensembl species code.
#' @param genome BSgenome object.
#' @param nuclease Nuclease type.
#' @param methods Scoring methods.
#' @param modern_methods Modern method names.
#' @param legacy_methods Legacy method names.
#' @param skip_plots Whether to skip plot generation.
#'
#' @return Named list of results (or error objects for failed genes).
#'
#' @noRd
process_batch_parallel <- function(gene_ids,
                                   species,
                                   genome,
                                   nuclease,
                                   methods,
                                   modern_methods,
                                   legacy_methods,
                                   skip_plots) {

  # Capture the Python path from the parent environment
  python_path <- Sys.getenv("RETICULATE_PYTHON")

  # Determine BSgenome package name for future.packages
  genome_pkg <- class(genome)[1]

  # Worker function
  process_single_gene_safe <- function(gene_id) {

    # ══════════════════════════════════════════════════════════════
    # Ensure correct Python environment
    # ══════════════════════════════════════════════════════════════

    if (!reticulate::py_available(initialize = FALSE)) {
      if (nzchar(python_path) && file.exists(python_path)) {
        Sys.setenv(RETICULATE_PYTHON = python_path)
      }
    }

    env_ready <- tryCatch({
      if (!mutateR::check_mutater_env()) {
        mutateR::activate_mutater_env()
      } else {
        TRUE
      }
    }, error = function(e) {
      FALSE
    })

    if (!env_ready) {
      if (reticulate::py_available()) {
        current_py <- tryCatch(reticulate::py_config()$python, error = function(e) "")
        if (!grepl("r-mutater", current_py)) {
          return(list(
            error = TRUE,
            message = paste0("Wrong Python environment: ", current_py),
            gene_id = gene_id
          ))
        }
      } else {
        return(list(
          error = TRUE,
          message = "Failed to initialize Python environment in worker",
          gene_id = gene_id
        ))
      }
    }

    # ══════════════════════════════════════════════════════════════
    # Process gene (with suppressed messages)
    # ══════════════════════════════════════════════════════════════

    tryCatch({
      # Suppress verbose messages from mutateR scoring functions
      result <- suppressMessages(suppressWarnings({
        .analyze_single_gene(
          gene_id = gene_id,
          species = species,
          genome = genome,
          nuclease = nuclease,
          methods = methods,
          modern_methods = modern_methods,
          legacy_methods = legacy_methods,
          transcript_id = NULL,
          id_type = "symbol",
          tracr = "Chen2013",
          deephf_var = "wt_u6",
          cor_method = "spearman",
          top_n = c(10, 20, 50),
          generate_plots = !skip_plots,
          quiet = TRUE
        )
      }))

      if (skip_plots) {
        result$plots <- NULL
      }

      return(result)

    }, error = function(e) {
      list(
        error = TRUE,
        message = e$message,
        gene_id = gene_id
      )
    })
  }

  # Run in parallel
  results <- future.apply::future_lapply(
    gene_ids,
    process_single_gene_safe,
    future.seed = TRUE,
    future.packages = c("mutateR", genome_pkg, "reticulate")
  )

  names(results) <- gene_ids
  return(results)
}


#' Run staged validation analysis
#'
#' Convenience function to run progressively larger analyses to validate
#' pipeline stability and estimate resource requirements before committing
#' to full exome-scale analysis.
#'
#' @param all_genes Character vector of all gene symbols to sample from.
#' @param species Ensembl species code.
#' @param genome BSgenome object.
#' @param base_output_dir Base directory for staged outputs.
#' @param stages Named integer vector of stage sizes (default: 20, 200, 2000).
#' @param n_workers Number of parallel workers.
#' @param seed Random seed for reproducible sampling.
#'
#' @return List of paths to results from each stage.
#'
#' @examples
#' \dontrun{
#' all_genes <- get_protein_coding_genes("hsapiens")
#'
#' staged_results <- run_staged_validation(
#'   all_genes = all_genes,
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   stages = c(quick = 20, small = 200, medium = 2000)
#' )
#' }
#'
#' @export
run_staged_validation <- function(all_genes,
                                  species,
                                  genome,
                                  base_output_dir = "staged_analysis",
                                  stages = c(quick = 20, small = 200, medium = 2000),
                                  n_workers = 8,
                                  seed = 42) {

  dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

  results <- list()

  for (stage_name in names(stages)) {
    stage_size <- stages[[stage_name]]

    message("\n")
    message("═══════════════════════════════════════════════════════════════")
    message("Stage: ", stage_name, " (", stage_size, " genes)")
    message("═══════════════════════════════════════════════════════════════\n")

    # Sample genes reproducibly
    set.seed(seed)
    stage_genes <- sample(all_genes, min(stage_size, length(all_genes)))

    # Run analysis
    output_dir <- file.path(base_output_dir, stage_name)

    result_path <- run_exome_analysis(
      gene_ids = stage_genes,
      species = species,
      genome = genome,
      output_dir = output_dir,
      n_workers = n_workers,
      batch_size = min(50, stage_size),
      skip_plots = TRUE
    )

    results[[stage_name]] <- result_path

    # Brief summary
    if (file.exists(result_path)) {
      stage_result <- readRDS(result_path)
      message("\nStage '", stage_name, "' complete: ",
              length(stage_result), " genes processed")
    }
  }

  return(results)
}
