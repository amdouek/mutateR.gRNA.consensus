#' Run consensus scoring analysis for one or more genes
#'
#' Wrapper function that executes the full consensus scoring exploration pipeline
#' for a single gene or a vector of genes. Handles the full workflow:
#' gene lookup → exon retrieval → gRNA finding → multi-method scoring → agreement analysis.
#'
#' @param gene_ids Character vector. One or more gene symbols (e.g., c("TP53", "BRCA1"))
#'        or Ensembl gene IDs.
#' @param species Character. Ensembl species code (e.g., "hsapiens", "mmusculus", "drerio").
#' @param genome BSgenome object for the target species.
#' @param nuclease Character. One of "Cas9" (default), "Cas12a", or "enCas12a".
#' @param methods Character vector. Scoring methods to use. If NULL (default), auto-selects:
#'        - Cas9: c("ruleset1", "deepspcas9", "ruleset3", "deephf")
#'        - Cas12a: c("deepcpf1", "enpamgb")
#' @param transcript_ids Named character vector. Optional specific transcript IDs keyed by gene_id.
#'        E.g., c("TP53" = "ENST00000269305", "BRCA1" = "ENST00000357654")
#' @param id_type Character. One of "symbol" (default) or "ensembl_gene_id".
#' @param tracr Character. tracrRNA for RuleSet3 (default "Chen2013").
#' @param deephf_var Character. DeepHF variant (default "wt_u6").
#' @param cor_method Character. Correlation method: "spearman" (default) or "pearson".
#' @param top_n Integer vector. Thresholds for overlap analysis (default: c(10, 20, 50)).
#' @param generate_plots Logical. Whether to generate plots (default TRUE).
#' @param quiet Logical. Suppress progress messages (default FALSE).
#' @param stop_on_error Logical. If FALSE (default), continues processing remaining genes
#'        when one fails. If TRUE, stops immediately on first error.
#'
#' @return
#'   - If single gene: A "mutateR_consensus_analysis" object
#'   - If multiple genes: A "mutateR_consensus_batch" object (named list of individual analyses)
#'     with additional batch-level methods and comparison utilities
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Single gene
#' tp53 <- analyze_score_consensus("TP53", "hsapiens", BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Multiple genes
#' results <- analyze_score_consensus(
#'   c("TP53", "BRCA1", "EGFR", "CD4", "MYC"),
#'   "hsapiens",
#'   BSgenome.Hsapiens.UCSC.hg38
#' )
#'
#' # Access individual results
#' results$TP53
#' results$BRCA1$plots$matrix
#'
#' # Get comparison table
#' summary(results)
#'
#' # Large-scale analysis (all protein-coding genes - use with caution!)
#' # all_genes <- get_all_protein_coding_genes("hsapiens")
#' # all_results <- analyze_score_consensus(all_genes, "hsapiens", genome, quiet = TRUE)
#' }
#'
#' @export
analyze_score_consensus <- function(gene_ids,
                                    species,
                                    genome,
                                    nuclease = c("Cas9", "Cas12a", "enCas12a"),
                                    methods = NULL,
                                    transcript_ids = NULL,
                                    id_type = c("symbol", "ensembl_gene_id"),
                                    tracr = "Chen2013",
                                    deephf_var = "wt_u6",
                                    cor_method = c("spearman", "pearson"),
                                    top_n = c(10, 20, 50),
                                    generate_plots = TRUE,
                                    quiet = FALSE,
                                    stop_on_error = FALSE) {

  # --- Argument matching ---
  nuclease <- match.arg(nuclease)
  id_type <- match.arg(id_type)
  cor_method <- match.arg(cor_method)

  if (!inherits(genome, "BSgenome")) {
    stop("Please supply a valid BSgenome object.")
  }

  # --- Set default methods based on nuclease ---
  if (is.null(methods)) {
    methods <- switch(nuclease,
                      "Cas9" = c("ruleset1", "deepspcas9", "ruleset3", "deephf"),
                      "Cas12a" = c("deepcpf1", "enpamgb"),
                      "enCas12a" = c("deepcpf1", "enpamgb")
    )
  }

  # --- Set method classification ---
  if (nuclease == "Cas9") {
    modern_methods <- intersect(c("deepspcas9", "deephf", "ruleset3"), methods)
    legacy_methods <- intersect(c("ruleset1"), methods)
  } else {
    modern_methods <- methods
    legacy_methods <- character(0)
  }

  # --- Single gene: delegate to internal function ---
  if (length(gene_ids) == 1) {
    tx_id <- if (!is.null(transcript_ids)) transcript_ids[[gene_ids]] else NULL
    return(.analyze_single_gene(
      gene_id = gene_ids,
      species = species,
      genome = genome,
      nuclease = nuclease,
      methods = methods,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      transcript_id = tx_id,
      id_type = id_type,
      tracr = tracr,
      deephf_var = deephf_var,
      cor_method = cor_method,
      top_n = top_n,
      generate_plots = generate_plots,
      quiet = quiet
    ))
  }

  # --- Multiple genes: batch processing ---
  n_genes <- length(gene_ids)
  if (!quiet) {
    message("╔═══════════════════════════════════════════════════════════════╗")
    message("║         mutateR Batch Consensus Analysis                      ║")
    message("╚═══════════════════════════════════════════════════════════════╝")
    message("\nProcessing ", n_genes, " genes for ", species, " (", nuclease, ")")
    message("Methods: ", paste(methods, collapse = ", "), "\n")
  }

  results <- list()
  failed_genes <- character()

  for (i in seq_along(gene_ids)) {
    gene <- gene_ids[i]

    if (!quiet) {
      message(sprintf("[%d/%d] Processing %s...", i, n_genes, gene))
    }

    tx_id <- if (!is.null(transcript_ids)) transcript_ids[[gene]] else NULL

    result <- tryCatch({
      .analyze_single_gene(
        gene_id = gene,
        species = species,
        genome = genome,
        nuclease = nuclease,
        methods = methods,
        modern_methods = modern_methods,
        legacy_methods = legacy_methods,
        transcript_id = tx_id,
        id_type = id_type,
        tracr = tracr,
        deephf_var = deephf_var,
        cor_method = cor_method,
        top_n = top_n,
        generate_plots = generate_plots,
        quiet = TRUE  # Suppress individual gene messages in batch mode
      )
    }, error = function(e) {
      if (stop_on_error) {
        stop("Error processing '", gene, "': ", e$message)
      }
      if (!quiet) message("Failed: ", e$message)
      failed_genes <<- c(failed_genes, gene)
      return(NULL)
    })

    if (!is.null(result)) {
      results[[gene]] <- result
    }
  }

  # --- Build batch result object ---
  batch_result <- structure(
    results,
    class = c("mutateR_consensus_batch", "list"),
    metadata = list(
      gene_ids = gene_ids,
      n_requested = n_genes,
      n_successful = length(results),
      n_failed = length(failed_genes),
      failed_genes = failed_genes,
      species = species,
      nuclease = nuclease,
      methods = methods,
      modern_methods = if (nuclease == "Cas9") {
        intersect(c("deepspcas9", "deephf", "ruleset3"), methods)
      } else methods,
      legacy_methods = if (nuclease == "Cas9") {
        intersect(c("ruleset1"), methods)
      } else character(0),
      cor_method = cor_method,
      timestamp = Sys.time()
    )
  )

  if (!quiet) {
    message("\n───────────────────────────────────────────────────────────────")
    message("Batch complete: ", length(results), "/", n_genes, " genes successful")
    if (length(failed_genes) > 0) {
      message("Failed genes: ", paste(failed_genes, collapse = ", "))
    }
    message("───────────────────────────────────────────────────────────────\n")
  }

  return(batch_result)
}


#' Internal function: Analyse a single gene
#' @noRd
.analyze_single_gene <- function(gene_id,
                                 species,
                                 genome,
                                 nuclease,
                                 methods,
                                 modern_methods,
                                 legacy_methods,
                                 transcript_id,
                                 id_type,
                                 tracr,
                                 deephf_var,
                                 cor_method,
                                 top_n,
                                 generate_plots,
                                 quiet) {

  # --- Initialise result structure ---
  result <- list(
    metadata = list(
      gene_id = gene_id,
      species = species,
      nuclease = nuclease,
      methods = methods,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      id_type = id_type,
      tracr = tracr,
      deephf_var = deephf_var,
      cor_method = cor_method,
      timestamp = Sys.time()
    ),
    exons = NULL,
    sites = NULL,
    scores = NULL,
    agreement = NULL,
    plots = NULL
  )

  # --- Step 1: Gene/Transcript Information ---
  if (!quiet) message("═══ Analysing ", gene_id, " (", species, ", ", nuclease, ") ═══")
  if (!quiet) message("\n[1/5] Retrieving gene information...")

  tx_info <- tryCatch({
    get_gene_info(gene_id, species, id_type = id_type)
  }, error = function(e) {
    stop("Failed to retrieve gene info for '", gene_id, "': ", e$message)
  })

  # Determine transcript to use
  if (is.null(transcript_id)) {
    if (!is.null(tx_info$canonical) && nrow(tx_info$canonical) > 0) {
      tx_id <- tx_info$canonical$ensembl_transcript_id[1]
      gene_symbol <- tx_info$canonical$external_gene_name[1]
    } else {
      tx_id <- tx_info$all$ensembl_transcript_id[1]
      gene_symbol <- tx_info$all$external_gene_name[1]
      if (!quiet) warning("No canonical transcript found; using first available: ", tx_id)
    }
  } else {
    tx_id <- transcript_id
    gene_symbol <- if (!is.null(tx_info$all)) {
      tx_info$all$external_gene_name[tx_info$all$ensembl_transcript_id == tx_id][1]
    } else {
      gene_id
    }
  }

  result$metadata$transcript_id <- tx_id
  result$metadata$gene_symbol <- gene_symbol

  if (!quiet) message("   Using transcript: ", tx_id, " (", gene_symbol, ")")

  # --- Step 2: Exon Structures ---
  if (!quiet) message("\n[2/5] Retrieving exon structures...")

  exons_gr <- tryCatch({
    get_exon_structures(tx_id, species, output = "GRanges")
  }, error = function(e) {
    stop("Failed to retrieve exon structures: ", e$message)
  })

  result$exons <- exons_gr
  n_exons <- length(exons_gr)

  if (!quiet) message("   Found ", n_exons, " exons")

  # --- Step 3: Find gRNA Sites ---
  if (!quiet) message("\n[3/5] Identifying ", nuclease, " target sites...")

  sites_gr <- tryCatch({
    switch(nuclease,
           "Cas9" = find_cas9_sites(exons_gr, genome, pam = "NGG", protospacer_length = 20),
           "Cas12a" = find_cas12a_sites(exons_gr, genome, pam = "TTTV", protospacer_length = 23),
           "enCas12a" = find_cas12a_sites(exons_gr, genome, pam = "TTTN", protospacer_length = 23)
    )
  }, error = function(e) {
    stop("Failed to find gRNA sites: ", e$message)
  })

  if (is.null(sites_gr) || length(sites_gr) == 0) {
    warning("No gRNA sites found for ", gene_id)
    result$metadata$n_sites <- 0
    class(result) <- c("mutateR_consensus_analysis", class(result))
    return(result)
  }

  result$sites <- sites_gr
  n_sites <- length(sites_gr)
  result$metadata$n_sites <- n_sites

  if (!quiet) message("   Found ", n_sites, " candidate sites")

  # --- Step 4: Multi-Method Scoring ---
  if (!quiet) message("\n[4/5] Scoring with ", length(methods), " methods...")

  multi_scores <- tryCatch({
    score_grnas_multi(
      sites_gr,
      methods = methods,
      tracr = tracr,
      deephf_var = deephf_var,
      quiet = quiet
    )
  }, error = function(e) {
    stop("Multi-method scoring failed: ", e$message)
  })

  result$scores <- multi_scores

  # --- Step 5: Agreement Analysis ---
  if (!quiet) message("\n[5/5] Computing agreement statistics...")

  # Check we have at least 2 methods with scores
  scored_methods <- methods[sapply(methods, function(m) {
    m %in% names(multi_scores) && sum(!is.na(multi_scores[[m]])) > 0
  })]

  if (length(scored_methods) < 2) {
    warning("Fewer than 2 methods returned valid scores. Skipping agreement analysis.")
    class(result) <- c("mutateR_consensus_analysis", class(result))
    return(result)
  }

  agreement <- tryCatch({
    suppressWarnings({
      summarize_score_agreement(
        multi_scores,
        methods = scored_methods,
        top_n = top_n,
        cor_method = cor_method
      )
    })
  }, error = function(e) {
    warning("Agreement analysis failed: ", e$message)
    NULL
  })

  result$agreement <- agreement

  # --- Step 6: Generate Plots (Optional) ---
  if (generate_plots && !is.null(agreement)) {
    if (!quiet) message("\n[+] Generating plots...")

    result$plots <- list()

    # Correlation matrix
    result$plots$matrix <- tryCatch({
      plot_score_correlations(multi_scores, methods = scored_methods,
                              plot_type = "matrix", cor_method = cor_method)
    }, error = function(e) {
      if (!quiet) message("   Matrix plot failed: ", e$message)
      NULL
    })

    # Scatter plots
    result$plots$scatter <- tryCatch({
      plot_score_correlations(multi_scores, methods = scored_methods,
                              plot_type = "scatter", cor_method = cor_method,
                              highlight_top = 20)
    }, error = function(e) {
      if (!quiet) message("   Scatter plot failed: ", e$message)
      NULL
    })
  }

  # --- Finalise ---
  class(result) <- c("mutateR_consensus_analysis", class(result))

  if (!quiet) {
    message("\n═══ Analysis complete ═══")
    message("Access results via: $metadata, $exons, $sites, $scores, $agreement, $plots")
  }

  return(result)
}


#' Print method for batch consensus analysis
#' @export
print.mutateR_consensus_batch <- function(x, ...) {

  meta <- attr(x, "metadata")

  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║         mutateR Batch Consensus Analysis                      ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

  cat("Species:    ", meta$species, "\n")
  cat("Nuclease:   ", meta$nuclease, "\n")
  cat("Methods:    ", paste(meta$methods, collapse = ", "), "\n")
  cat("Correlation:", meta$cor_method, "\n\n")

  cat("Genes requested:  ", meta$n_requested, "\n")
  cat("Genes successful: ", meta$n_successful, "\n")
  cat("Genes failed:     ", meta$n_failed, "\n")

  if (meta$n_failed > 0) {
    cat("Failed genes:     ", paste(meta$failed_genes, collapse = ", "), "\n")
  }

  cat("\n")
  cat("───────────────────────────────────────────────────────────────\n")
  cat("Use summary(x) for comparison table across genes.\n")
  cat("Access individual results via x$<gene_name> or x[[\"gene_name\"]].\n")
  cat("───────────────────────────────────────────────────────────────\n\n")

  invisible(x)
}


#' Summary method for batch consensus analysis
#' @export
summary.mutateR_consensus_batch <- function(object, ...) {
  compare_consensus_analyses(object)
}


#' Extract comparison table from batch analysis
#'
#' @param x A mutateR_consensus_batch object
#' @param include_correlations Logical. Include individual pairwise correlations (default FALSE).
#'
#' @return A data frame with one row per gene and summary statistics.
#' @export
as.data.frame.mutateR_consensus_batch <- function(x,
                                                  row.names = NULL,
                                                  optional = FALSE,
                                                  include_correlations = FALSE,
                                                  ...) {

  if (length(x) == 0) {
    return(data.frame())
  }

  meta <- attr(x, "metadata")

  rows <- lapply(names(x), function(gene) {
    obj <- x[[gene]]

    row_data <- data.frame(
      gene = gene,
      gene_symbol = obj$metadata$gene_symbol,
      transcript_id = obj$metadata$transcript_id,
      n_exons = if (!is.null(obj$exons)) length(obj$exons) else NA_integer_,
      n_sites = obj$metadata$n_sites,
      stringsAsFactors = FALSE
    )

    if (!is.null(obj$agreement)) {
      cor_df <- obj$agreement$correlations

      row_data$mean_correlation <- mean(cor_df$correlation, na.rm = TRUE)
      row_data$min_correlation <- min(cor_df$correlation, na.rm = TRUE)
      row_data$max_correlation <- max(cor_df$correlation, na.rm = TRUE)

      # Individual correlations (optional)
      if (include_correlations) {
        for (i in seq_len(nrow(cor_df))) {
          col_name <- paste0("cor_", cor_df$method1[i], "_", cor_df$method2[i])
          row_data[[col_name]] <- cor_df$correlation[i]
        }
      }

      # Enrichment statistics
      for (n in unique(obj$agreement$rank_overlaps$top_n)) {
        enrich <- obj$agreement$rank_overlaps[obj$agreement$rank_overlaps$top_n == n, ]
        row_data[[paste0("enrich_top", n)]] <- mean(enrich$enrichment, na.rm = TRUE)
      }

      # Concordance metrics
      if (!is.null(obj$agreement$concordant_top) && nrow(obj$agreement$concordant_top) > 0) {
        row_data$best_max_rank <- min(obj$agreement$concordant_top$max_rank, na.rm = TRUE)
        row_data$n_top50_concordant <- sum(obj$agreement$concordant_top$max_rank <= 50, na.rm = TRUE)
      } else {
        row_data$best_max_rank <- NA_integer_
        row_data$n_top50_concordant <- 0L
      }

    } else {
      row_data$mean_correlation <- NA_real_
      row_data$min_correlation <- NA_real_
      row_data$max_correlation <- NA_real_
      row_data$best_max_rank <- NA_integer_
      row_data$n_top50_concordant <- NA_integer_
    }

    return(row_data)
  })

  # Combine with consistent columns
  all_cols <- unique(unlist(lapply(rows, names)))
  result <- do.call(rbind, lapply(rows, function(row) {
    missing <- setdiff(all_cols, names(row))
    for (col in missing) row[[col]] <- NA
    row[, all_cols, drop = FALSE]
  }))

  rownames(result) <- NULL
  return(result)
}
