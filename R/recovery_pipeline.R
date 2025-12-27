#' Check if seqlevels are on primary assembly
#'
#' Determines whether a GRanges object contains only standard chromosomal
#' sequences (primary assembly) or includes patches/alternative haplotypes.
#'
#' @param gr GRanges object
#' @param style Character. Expected seqlevel style: "Ensembl" or "UCSC"
#'
#' @return List with components:
#'   \describe{
#'     \item{is_primary}{Logical. TRUE if all seqlevels are primary assembly}
#'     \item{seqlevels_used}{Character vector of seqlevels containing ranges}
#'     \item{nonstandard}{Character vector of non-primary seqlevels (empty if none)}
#'     \item{categories}{Named character vector categorising each non-standard seqlevel}
#'   }
#'
#' @details
#' Seqlevel categories:
#' \describe{
#'   \item{fix_patch}{GRC fix patch (corrects primary assembly error)}
#'   \item{novel_patch}{GRC novel patch (new sequence not in primary)}
#'   \item{mhc_haplotype}{MHC region alternative haplotype}
#'   \item{alt_haplotype}{Other alternative haplotype}
#'   \item{unknown}{Unrecognised non-standard sequence}
#' }
#'
#' @noRd
check_seqlevel_status <- function(gr, style = c("Ensembl", "UCSC")){

  style <- match.arg(style)

  # Define standard chromosomes
  if (style == "Ensembl") {
    standard_seqlevels <- c(as.character(1:22), "X", "Y", "MT")
  } else {
    standard_seqlevels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  }

  # Get seqlevels actually in use
  seqlevels_used <- unique(as.character(GenomicRanges::seqnames(gr)))

  # Identify non-standard
  nonstandard <- setdiff(seqlevels_used, standard_seqlevels)

  # Categorise non-standard seqlevels
  categories <- character()
  if (length(nonstandard) > 0) {
    categories <- vapply(nonstandard, function(sl) {
      sl_upper <- toupper(sl)
      if (grepl("_MHC_", sl_upper)) {
        return("mhc_haplotype")
      } else if (grepl("_PATCH$", sl_upper)) {
        # Could potentially distinguish FIX vs NOVEL by querying Ensembl
        # For now, categorise as generic patch
        return("patch")
      } else if (grepl("^HSCHR", sl_upper) || grepl("_CTG", sl_upper)) {
        return("alt_haplotype")
      } else if (grepl("^CHR_", sl_upper)) {
        return("alt_haplotype")
      } else {
        return("unknown")
      }
    }, character(1))
    names(categories) <- nonstandard
  }

  list(
    is_primary = length(nonstandard) == 0,
    seqlevels_used = seqlevels_used,
    nonstandard = nonstandard,
    categories = categories
  )
}

#' Generate structured failure message for non-standard seqlevels
#'
#' Creates a parseable error message containing seqlevel information
#' for downstream recovery processing.
#'
#' @param seqlevel_status Output from check_seqlevel_status()
#' @return Character string with structured failure information
#'
#' @noRd
format_seqlevel_failure <- function(seqlevel_status) {

  if (seqlevel_status$is_primary) {
    return(NULL)
  }

  nonstandard <- seqlevel_status$nonstandard
  categories <- seqlevel_status$categories

  # Format: "Non-standard seqlevel: <seqlevel> (<category>)"
  # Use first non-standard seqlevel for primary message
  primary_seqlevel <- nonstandard[1]
  primary_category <- categories[primary_seqlevel]

  msg <- sprintf(
    "Non-standard seqlevel: %s (%s)",
    primary_seqlevel,
    primary_category
  )

  # Add additional seqlevels if present

  if (length(nonstandard) > 1) {
    msg <- paste0(msg, sprintf(" [+%d others]", length(nonstandard) - 1))
  }

  return(msg)
}

#' Retrieve sequence from Ensembl REST API
#'
#' Fetches DNA sequence for a specified genomic region using the Ensembl
#' REST API. Works for any annotated sequence including patches and
#' alternative haplotypes.
#'
#' @param seqlevel Character. Sequence name (e.g., "1", "HG401_PATCH")
#' @param start Integer. Start coordinate (1-based)
#' @param end Integer. End coordinate (1-based, inclusive)
#' @param strand Character. "+" or "-"
#' @param species Character. Species name for API (default "human")
#' @param expand Integer. Bases to expand on each side (default 0)
#'
#' @return Character string containing the DNA sequence
#'
#' @details
#' Uses the Ensembl REST API endpoint:
#' \code{GET /sequence/region/:species/:region}
#'
#' Rate limiting: Ensembl recommends maximum 15 requests/second.
#' This function includes a small delay to stay within limits when
#' called repeatedly.
#'
#' @noRd
fetch_ensembl_sequence <- function(seqlevel,
                                   start,
                                   end,
                                   strand = "+",
                                   species = "human",
                                   expand = 0) {

  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' required for Ensembl API calls. ",
         "Install with: install.packages('httr')")
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' required for Ensembl API calls. ",
         "Install with: install.packages('jsonlite')")
  }

  # Convert strand to Ensembl format
  strand_num <- if (strand == "-") -1 else 1

  # Apply expansion
  query_start <- max(1, start - expand)
  query_end <- end + expand

  # Build region string
  region <- sprintf("%s:%d..%d:%d", seqlevel, query_start, query_end, strand_num)

  # Build URL
  base_url <- "https://rest.ensembl.org"
  endpoint <- sprintf("/sequence/region/%s/%s", species, region)
  url <- paste0(base_url, endpoint)

  # Make request
  response <- tryCatch({
    httr::GET(
      url,
      httr::content_type("text/plain"),
      httr::accept("text/plain"),
      httr::timeout(30)
    )
  }, error = function(e) {
    stop("Ensembl API request failed: ", e$message)
  })

  # Check status
  if (httr::status_code(response) != 200) {
    # Try to get error message
    error_content <- tryCatch({
      httr::content(response, as = "text", encoding = "UTF-8")
    }, error = function(e) "Unknown error")

    stop(sprintf("Ensembl API error (HTTP %d) for region %s: %s",
                 httr::status_code(response), region, error_content))
  }

  # Extract sequence
  sequence <- httr::content(response, as = "text", encoding = "UTF-8")
  sequence <- gsub("\\s", "", sequence)
  sequence <- toupper(sequence)

  # Small delay to respect rate limiting
  Sys.sleep(0.07)  # ~14 requests/second max

  return(sequence)
}


#' Retrieve sequences for multiple exons
#'
#' Batch retrieval of exon sequences from Ensembl REST API with
#' appropriate context windows for gRNA finding.
#'
#' @param exons_gr GRanges object containing exon coordinates
#' @param species Character. Species for API (default "human")
#' @param context_upstream Integer. Bases upstream to include (default 4)
#' @param context_downstream Integer. Bases downstream to include (default 6)
#'
#' @return DNAStringSet with sequences, named by range index
#'
#' @noRd
fetch_exon_sequences <- function(exons_gr,
                                 species = "human",
                                 context_upstream = 4,
                                 context_downstream = 6) {

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' required.")
  }

  n_exons <- length(exons_gr)
  sequences <- character(n_exons)

  for (i in seq_len(n_exons)) {
    seqlevel <- as.character(GenomicRanges::seqnames(exons_gr)[i])
    start_pos <- GenomicRanges::start(exons_gr)[i]
    end_pos <- GenomicRanges::end(exons_gr)[i]
    strand_char <- as.character(GenomicRanges::strand(exons_gr)[i])

    # Fetch with context
    # Note: for negative strand, Ensembl returns reverse complement
    sequences[i] <- fetch_ensembl_sequence(
      seqlevel = seqlevel,
      start = start_pos,
      end = end_pos,
      strand = strand_char,
      species = species,
      expand = max(context_upstream, context_downstream)
    )
  }

  # Convert to DNAStringSet
  seq_set <- Biostrings::DNAStringSet(sequences)
  names(seq_set) <- paste0("exon_", seq_len(n_exons))

  return(seq_set)
}

#' Find Cas9 PAM sites in a sequence
#'
#' Scans a DNA sequence for NGG PAM motifs and extracts protospacer
#' sequences and context windows. Replicates the output format of
#' mutateR::find_cas9_sites() but operates on sequence strings rather
#' than BSgenome objects.
#'
#' @param sequence Character. DNA sequence to scan
#' @param seqlevel Character. Sequence name for coordinate reporting
#' @param seq_start Integer. Genomic start coordinate of the sequence
#' @param seq_strand Character. Strand of the sequence ("+" or "-")
#' @param exon_rank Integer. Exon number for annotation
#' @param pam Character. PAM sequence (default "NGG")
#' @param protospacer_length Integer. Length of protospacer (default 20)
#'
#' @return Data frame with columns matching find_cas9_sites() output:
#'   protospacer_sequence, sequence_context, chr, start, end, strand, exon_rank
#'
#' @noRd
scan_cas9_sites <- function(sequence,
                            seqlevel,
                            seq_start,
                            seq_strand,
                            exon_rank,
                            pam = "NGG",
                            protospacer_length = 20L) {

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' required.")
  }

  sequence <- toupper(sequence)
  seq_length <- nchar(sequence)

  # Context window: 4bp upstream + 20bp protospacer + 3bp PAM + 3bp downstream = 30bp
  context_upstream <- 4L
  context_downstream <- 3L
  context_length <- context_upstream + protospacer_length + 3L + context_downstream

  sites <- list()

  # --- Forward strand sites (NGG PAM) ---
  # PAM is at positions i, i+1, i+2 where sequence[i+1:i+2] == "GG"
  # Protospacer is 20nt upstream of PAM
  ngg_pattern <- "(?=(.GG))"
  forward_matches <- gregexpr(ngg_pattern, sequence, perl = TRUE)[[1]]

  if (forward_matches[1] != -1) {
    for (match_pos in forward_matches) {
      pam_start <- match_pos  # Position of N in NGG
      protospacer_start <- pam_start - protospacer_length
      protospacer_end <- pam_start - 1
      context_start <- protospacer_start - context_upstream
      context_end <- pam_start + 2 + context_downstream  # PAM is 3bp

      # Check bounds
      if (protospacer_start < 1) next
      if (context_end > seq_length) next
      if (context_start < 1) next

      protospacer <- substr(sequence, protospacer_start, protospacer_end)
      context <- substr(sequence, context_start, context_end)

      # Skip if protospacer contains N
      if (grepl("N", protospacer)) next

      # Calculate genomic coordinates
      if (seq_strand == "+") {
        genomic_start <- seq_start + protospacer_start - 1
        genomic_end <- seq_start + protospacer_end - 1
        site_strand <- "+"
      } else {
        # For negative strand input, coordinates are already flipped
        genomic_end <- seq_start + (seq_length - protospacer_start)
        genomic_start <- seq_start + (seq_length - protospacer_end)
        site_strand <- "+"
      }

      sites[[length(sites) + 1]] <- data.frame(
        protospacer_sequence = protospacer,
        sequence_context = context,
        chr = seqlevel,
        start = genomic_start,
        end = genomic_end,
        strand = site_strand,
        exon_rank = exon_rank,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Reverse strand sites (CCN on forward = NGG on reverse) ---
  ccn_pattern <- "(?=(CC.))"
  reverse_matches <- gregexpr(ccn_pattern, sequence, perl = TRUE)[[1]]

  if (reverse_matches[1] != -1) {
    for (match_pos in reverse_matches) {
      # PAM (CCN) is at match_pos:match_pos+2
      # Protospacer is 20nt downstream (on forward strand), which is upstream on reverse
      pam_end <- match_pos + 2
      protospacer_start <- pam_end + 1
      protospacer_end <- pam_end + protospacer_length
      context_start <- match_pos - context_downstream
      context_end <- protospacer_end + context_upstream

      # Check bounds
      if (protospacer_end > seq_length) next
      if (context_end > seq_length) next
      if (context_start < 1) next

      # Get sequences and reverse complement
      protospacer_fwd <- substr(sequence, protospacer_start, protospacer_end)
      context_fwd <- substr(sequence, context_start, context_end)

      protospacer <- as.character(Biostrings::reverseComplement(
        Biostrings::DNAString(protospacer_fwd)
      ))
      context <- as.character(Biostrings::reverseComplement(
        Biostrings::DNAString(context_fwd)
      ))

      # Skip if protospacer contains N
      if (grepl("N", protospacer)) next

      # Calculate genomic coordinates
      if (seq_strand == "+") {
        genomic_start <- seq_start + protospacer_start - 1
        genomic_end <- seq_start + protospacer_end - 1
        site_strand <- "-"
      } else {
        genomic_end <- seq_start + (seq_length - protospacer_start)
        genomic_start <- seq_start + (seq_length - protospacer_end)
        site_strand <- "-"
      }

      sites[[length(sites) + 1]] <- data.frame(
        protospacer_sequence = protospacer,
        sequence_context = context,
        chr = seqlevel,
        start = genomic_start,
        end = genomic_end,
        strand = site_strand,
        exon_rank = exon_rank,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(sites) == 0) {
    return(data.frame(
      protospacer_sequence = character(),
      sequence_context = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      strand = character(),
      exon_rank = integer(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, sites)
}


#' Find all Cas9 sites in exons using Ensembl sequences
#'
#' Main function to identify gRNA target sites from exon coordinates
#' by retrieving sequences from Ensembl REST API.
#'
#' @param exons_gr GRanges object containing exon coordinates
#' @param species Character. Species for Ensembl API (default "human")
#' @param pam Character. PAM sequence (default "NGG")
#' @param protospacer_length Integer. Protospacer length (default 20)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return GRanges object matching format of mutateR::find_cas9_sites()
#'
#' @noRd
find_sites_from_ensembl <- function(exons_gr,
                                    species = "human",
                                    pam = "NGG",
                                    protospacer_length = 20L,
                                    quiet = FALSE) {

  n_exons <- length(exons_gr)
  if (n_exons == 0) {
    stop("No exons provided")
  }

  if (!quiet) {
    message("Retrieving sequences from Ensembl for ", n_exons, " exons...")
  }

  # Context needed: 4bp upstream + 3bp downstream of protospacer+PAM region
  context_expand <- 10  # Conservative expansion

  all_sites <- list()

  for (i in seq_len(n_exons)) {
    seqlevel <- as.character(GenomicRanges::seqnames(exons_gr)[i])
    start_pos <- GenomicRanges::start(exons_gr)[i]
    end_pos <- GenomicRanges::end(exons_gr)[i]
    strand_char <- as.character(GenomicRanges::strand(exons_gr)[i])

    # Get exon rank from metadata if available
    exon_rank <- if ("exon_rank" %in% names(GenomicRanges::mcols(exons_gr))) {
      GenomicRanges::mcols(exons_gr)$exon_rank[i]
    } else {
      i
    }

    # Fetch sequence with context
    sequence <- tryCatch({
      fetch_ensembl_sequence(
        seqlevel = seqlevel,
        start = start_pos,
        end = end_pos,
        strand = "+",  # Always fetch forward, handle strand in scanning
        species = species,
        expand = context_expand
      )
    }, error = function(e) {
      if (!quiet) message("  Failed to fetch exon ", i, ": ", e$message)
      return(NULL)
    })

    if (is.null(sequence)) next

    # Scan for PAM sites
    # Adjust start position for context expansion
    adjusted_start <- start_pos - context_expand

    sites_df <- scan_cas9_sites(
      sequence = sequence,
      seqlevel = seqlevel,
      seq_start = adjusted_start,
      seq_strand = strand_char,
      exon_rank = exon_rank,
      pam = pam,
      protospacer_length = protospacer_length
    )

    if (nrow(sites_df) > 0) {
      all_sites[[length(all_sites) + 1]] <- sites_df
    }
  }

  if (length(all_sites) == 0) {
    if (!quiet) message("No gRNA sites found")
    return(GenomicRanges::GRanges())
  }

  # Combine all sites
  sites_df <- do.call(rbind, all_sites)

  # Remove duplicates (sites spanning exon boundaries might be found twice)
  sites_df <- unique(sites_df)

  if (!quiet) {
    message("Found ", nrow(sites_df), " candidate sites")
  }

  # Convert to GRanges
  sites_gr <- GenomicRanges::GRanges(
    seqnames = sites_df$chr,
    ranges = IRanges::IRanges(start = sites_df$start, end = sites_df$end),
    strand = sites_df$strand,
    protospacer_sequence = sites_df$protospacer_sequence,
    sequence_context = sites_df$sequence_context,
    exon_rank = sites_df$exon_rank
  )

  return(sites_gr)
}

#' Recover failed genes using Ensembl REST API
#'
#' Processes genes that failed primary pipeline analysis due to non-standard
#' sequence locations by retrieving sequences directly from Ensembl.
#'
#' @param failed_genes Character vector of gene symbols, or path to failed_genes.rds
#' @param species Character. Ensembl species code (default "hsapiens")
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a" (default "Cas9")
#' @param methods Character vector. Scoring methods to use (default: auto-select)
#' @param tracr Character. tracrRNA for RuleSet3 (default "Chen2013")
#' @param deephf_var Character. DeepHF variant (default "wt_u6")
#' @param exclude_mhc Logical. Exclude MHC-region genes (default FALSE)
#' @param quiet Logical. Suppress progress messages (default FALSE)
#'
#' @return List with components:
#'   \describe{
#'     \item{recovered}{mutateR_consensus_batch of successfully recovered genes}
#'     \item{still_failed}{Data frame of genes that could not be recovered}
#'     \item{excluded_mhc}{Character vector of excluded MHC genes (if exclude_mhc=TRUE)}
#'     \item{summary}{Recovery statistics}
#'   }
#'
#' @details
#' This function provides an alternative pathway for genes whose canonical
#' transcripts are annotated on GRC patches or alternative haplotypes.
#' Instead of using a BSgenome object, it retrieves sequences directly
#' from the Ensembl REST API, which has access to all annotated sequences.
#'
#' The function:
#' \enumerate{
#'   \item Re-retrieves gene/transcript information from Ensembl
#'   \item Fetches exon sequences via REST API
#'   \item Scans for PAM sites and extracts protospacers
#'   \item Scores using the standard multi-method pipeline
#'   \item Computes agreement statistics
#' }
#'
#' Genes on MHC alternative haplotypes can optionally be excluded, as gRNAs
#' designed against these sequences will only work for individuals carrying
#' that specific HLA haplotype.
#'
#' @section Rate Limiting:
#' The Ensembl REST API has rate limits (~15 requests/second). This function
#' includes appropriate delays to stay within limits. For large numbers of
#' failed genes, expect approximately 1-2 seconds per gene.
#'
#' @examples
#' \dontrun{
#' # After running primary analysis
#' primary_results <- run_exome_analysis(genes, "hsapiens", genome)
#'
#' # Recover failed genes
#' recovery <- recover_failed_genes(
#'   "exome_analysis/failed_genes.rds",
#'   species = "hsapiens"
#' )
#'
#' # View recovery statistics
#' recovery$summary
#'
#' # Merge with primary results
#' merged <- merge_batch_results(
#'   readRDS("exome_analysis/final_results.rds"),
#'   recovery$recovered
#' )
#' }
#'
#' @export
recover_failed_genes <- function(failed_genes,
                                 species = "hsapiens",
                                 nuclease = c("Cas9", "Cas12a", "enCas12a"),
                                 methods = NULL,
                                 tracr = "Chen2013",
                                 deephf_var = "wt_u6",
                                 exclude_mhc = FALSE,
                                 quiet = FALSE) {

  nuclease <- match.arg(nuclease)

  # Load failed genes if path provided
  if (length(failed_genes) == 1 && file.exists(failed_genes)) {
    failed_list <- readRDS(failed_genes)
    gene_symbols <- names(failed_list)

    if (!quiet) {
      message("Loaded ", length(gene_symbols), " failed genes from: ", failed_genes)
    }
  } else {
    gene_symbols <- as.character(failed_genes)
  }

  if (length(gene_symbols) == 0) {
    message("No failed genes to recover")
    return(list(
      recovered = structure(list(), class = c("mutateR_consensus_batch", "list")),
      still_failed = data.frame(),
      excluded_mhc = character(),
      summary = list(n_input = 0, n_recovered = 0, n_failed = 0, n_mhc = 0)
    ))
  }

  # Set default methods
  if (is.null(methods)) {
    methods <- switch(nuclease,
                      "Cas9" = c("ruleset1", "deepspcas9", "ruleset3", "deephf"),
                      "Cas12a" = c("deepcpf1", "enpamgb"),
                      "enCas12a" = c("deepcpf1", "enpamgb")
    )
  }

  # Method classification
  if (nuclease == "Cas9") {
    modern_methods <- intersect(c("deepspcas9", "deephf", "ruleset3"), methods)
    legacy_methods <- intersect(c("ruleset1"), methods)
  } else {
    modern_methods <- methods
    legacy_methods <- character(0)
  }

  # Map species to Ensembl REST API format
  species_map <- c(
    "hsapiens" = "human",
    "mmusculus" = "mouse",
    "drerio" = "zebrafish",
    "rnorvegicus" = "rat"
  )
  api_species <- species_map[species]
  if (is.na(api_species)) api_species <- species

  if (!quiet) {
    cat("\n")
    cat("\u2554", paste(rep("\u2550", 63), collapse = ""), "\u2557\n", sep = "")
    cat("\u2551         Recovery Pipeline for Failed Genes                    \u2551\n")
    cat("\u255A", paste(rep("\u2550", 63), collapse = ""), "\u255D\n\n", sep = "")
    cat("Genes to recover: ", length(gene_symbols), "\n")
    cat("Species:          ", species, " (", api_species, ")\n", sep = "")
    cat("Nuclease:         ", nuclease, "\n")
    cat("Methods:          ", paste(methods, collapse = ", "), "\n")
    cat("Exclude MHC:      ", exclude_mhc, "\n\n")
  }

  # Track results
  recovered_results <- list()
  still_failed <- list()
  excluded_mhc <- character()

  for (i in seq_along(gene_symbols)) {
    gene <- gene_symbols[i]

    if (!quiet) {
      message(sprintf("[%d/%d] Processing %s...", i, length(gene_symbols), gene))
    }

    result <- tryCatch({
      recover_single_gene(
        gene_symbol = gene,
        species = species,
        api_species = api_species,
        nuclease = nuclease,
        methods = methods,
        modern_methods = modern_methods,
        legacy_methods = legacy_methods,
        tracr = tracr,
        deephf_var = deephf_var,
        exclude_mhc = exclude_mhc,
        quiet = TRUE
      )
    }, error = function(e) {
      list(success = FALSE, error = e$message, mhc_excluded = FALSE)
    })

    if (isTRUE(result$mhc_excluded)) {
      excluded_mhc <- c(excluded_mhc, gene)
      if (!quiet) message("  Excluded (MHC region)")
    } else if (isTRUE(result$success)) {
      recovered_results[[gene]] <- result$data
      if (!quiet) message("  Recovered: ", result$n_sites, " gRNA sites")
    } else {
      still_failed[[gene]] <- result$error
      if (!quiet) message("  Failed: ", result$error)
    }
  }

  # Build batch object for recovered genes
  recovered_batch <- structure(
    recovered_results,
    class = c("mutateR_consensus_batch", "list"),
    metadata = list(
      gene_ids = names(recovered_results),
      n_successful = length(recovered_results),
      n_failed = length(still_failed),
      failed_genes = names(still_failed),
      species = species,
      nuclease = nuclease,
      methods = methods,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      cor_method = "spearman",
      recovery_mode = TRUE,
      timestamp = Sys.time()
    )
  )

  # Build still_failed data frame
  still_failed_df <- if (length(still_failed) > 0) {
    data.frame(
      gene = names(still_failed),
      error = unlist(still_failed),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(gene = character(), error = character(), stringsAsFactors = FALSE)
  }

  # Summary
  summary_list <- list(
    n_input = length(gene_symbols),
    n_recovered = length(recovered_results),
    n_failed = length(still_failed),
    n_mhc_excluded = length(excluded_mhc),
    recovery_rate = length(recovered_results) / length(gene_symbols)
  )

  if (!quiet) {
    cat("\n")
    cat(paste(rep("\u2500", 65), collapse = ""), "\n")
    cat("Recovery Summary\n")
    cat(paste(rep("\u2500", 65), collapse = ""), "\n")
    cat("Input genes:      ", summary_list$n_input, "\n")
    cat("Recovered:        ", summary_list$n_recovered, "\n")
    cat("Still failed:     ", summary_list$n_failed, "\n")
    cat("MHC excluded:     ", summary_list$n_mhc_excluded, "\n")
    cat("Recovery rate:    ", sprintf("%.1f%%", 100 * summary_list$recovery_rate), "\n")
    cat(paste(rep("\u2500", 65), collapse = ""), "\n\n")
  }

  return(list(
    recovered = recovered_batch,
    still_failed = still_failed_df,
    excluded_mhc = excluded_mhc,
    summary = summary_list
  ))
}


#' Recover a single failed gene
#'
#' Internal function that performs the recovery process for one gene.
#'
#' @noRd
recover_single_gene <- function(gene_symbol,
                                species,
                                api_species,
                                nuclease,
                                methods,
                                modern_methods,
                                legacy_methods,
                                tracr,
                                deephf_var,
                                exclude_mhc,
                                quiet) {

  # Step 1: Get gene/transcript info
  tx_info <- tryCatch({
    mutateR::get_gene_info(gene_symbol, species, id_type = "symbol")
  }, error = function(e) {
    stop("Gene info retrieval failed: ", e$message)
  })

  # Select transcript
  if (!is.null(tx_info$canonical) && nrow(tx_info$canonical) > 0) {
    tx_id <- tx_info$canonical$ensembl_transcript_id[1]
    gene_sym <- tx_info$canonical$external_gene_name[1]
  } else {
    tx_id <- tx_info$all$ensembl_transcript_id[1]
    gene_sym <- tx_info$all$external_gene_name[1]
  }

  # Step 2: Get exon structures
  exons_gr <- tryCatch({
    mutateR::get_exon_structures(tx_id, species, output = "GRanges")
  }, error = function(e) {
    stop("Exon retrieval failed: ", e$message)
  })

  if (length(exons_gr) == 0) {
    stop("No exons found")
  }

  # Step 3: Check seqlevel status
  seqlevel_status <- check_seqlevel_status(exons_gr, style = "Ensembl")

  # Check for MHC exclusion
  if (exclude_mhc && any(seqlevel_status$categories == "mhc_haplotype")) {
    return(list(success = FALSE, mhc_excluded = TRUE))
  }

  # Step 4: Find gRNA sites via Ensembl
  sites_gr <- find_sites_from_ensembl(
    exons_gr = exons_gr,
    species = api_species,
    pam = if (nuclease == "Cas9") "NGG" else "TTTV",
    protospacer_length = if (nuclease == "Cas9") 20L else 23L,
    quiet = quiet
  )

  if (length(sites_gr) == 0) {
    stop("No gRNA sites found")
  }

  # Step 5: Score with multi-method pipeline
  multi_scores <- tryCatch({
    suppressMessages(suppressWarnings({
    score_grnas_multi(
      sites_gr,
      methods = methods,
      tracr = tracr,
      deephf_var = deephf_var,
      quiet = TRUE
    )
  }))
    }, error = function(e) {
    stop("Scoring failed: ", e$message)
  })

  # Step 6: Compute agreement
  agreement <- tryCatch({
    suppressWarnings({
      summarize_score_agreement(
        multi_scores,
        methods = methods,
        modern_methods = modern_methods,
        legacy_methods = legacy_methods,
        cor_method = "spearman"
      )
    })
  }, error = function(e) {
    NULL
  })

  # Step 7: Build result object
  result <- list(
    metadata = list(
      gene_id = gene_symbol,
      gene_symbol = gene_sym,
      transcript_id = tx_id,
      species = species,
      nuclease = nuclease,
      methods = methods,
      modern_methods = modern_methods,
      legacy_methods = legacy_methods,
      n_sites = length(sites_gr),
      seqlevel = seqlevel_status$seqlevels_used,
      seqlevel_category = if (length(seqlevel_status$categories) > 0) {
        seqlevel_status$categories[1]
      } else {
        "primary"
      },
      recovery_mode = TRUE,
      timestamp = Sys.time()
    ),
    exons = exons_gr,
    sites = sites_gr,
    scores = multi_scores,
    agreement = agreement,
    plots = NULL  # Skip plots in recovery mode
  )

  class(result) <- c("mutateR_consensus_analysis", class(result))

  return(list(
    success = TRUE,
    data = result,
    n_sites = length(sites_gr),
    mhc_excluded = FALSE
  ))
}

#' Merge primary and recovered batch results
#'
#' Combines results from the primary analysis pipeline with recovered
#' genes into a single mutateR_consensus_batch object.
#'
#' @param primary_batch Primary analysis results (mutateR_consensus_batch)
#' @param recovered_batch Recovered genes (from recover_failed_genes()$recovered)
#' @param mark_recovered Logical. Add metadata flag for recovered genes (default TRUE)
#'
#' @return Combined mutateR_consensus_batch object
#'
#' @examples
#' \dontrun{
#' primary <- readRDS("exome_analysis/final_results.rds")
#' recovery <- recover_failed_genes("exome_analysis/failed_genes.rds")
#'
#' merged <- merge_batch_results(primary, recovery$recovered)
#' }
#'
#' @export
merge_batch_results <- function(primary_batch,
                                recovered_batch,
                                mark_recovered = TRUE) {

  # Validate inputs
  if (!inherits(primary_batch, "mutateR_consensus_batch")) {
    stop("primary_batch must be a mutateR_consensus_batch object")
  }

  if (!inherits(recovered_batch, "mutateR_consensus_batch")) {
    stop("recovered_batch must be a mutateR_consensus_batch object")
  }

  primary_meta <- attr(primary_batch, "metadata")
  recovered_meta <- attr(recovered_batch, "metadata")

  # Check for overlapping genes (shouldn't happen)
  overlap <- intersect(names(primary_batch), names(recovered_batch))
  if (length(overlap) > 0) {
    warning("Overlapping genes found and will use primary version: ",
            paste(overlap, collapse = ", "))
    recovered_batch <- recovered_batch[!names(recovered_batch) %in% overlap]
  }

  # Add recovery flag if requested
  if (mark_recovered) {
    for (gene in names(recovered_batch)) {
      recovered_batch[[gene]]$metadata$was_recovered <- TRUE
    }
  }

  # Combine
  merged <- c(as.list(primary_batch), as.list(recovered_batch))

  # Build merged metadata
  merged_meta <- list(
    gene_ids = c(primary_meta$gene_ids, recovered_meta$gene_ids),
    n_successful = primary_meta$n_successful + recovered_meta$n_successful,
    n_failed = recovered_meta$n_failed,  # Only remaining failures
    n_recovered = recovered_meta$n_successful,
    failed_genes = recovered_meta$failed_genes,
    species = primary_meta$species,
    nuclease = primary_meta$nuclease,
    methods = primary_meta$methods,
    modern_methods = primary_meta$modern_methods,
    legacy_methods = primary_meta$legacy_methods,
    cor_method = primary_meta$cor_method,
    primary_timestamp = primary_meta$timestamp,
    recovery_timestamp = recovered_meta$timestamp,
    merged_timestamp = Sys.time()
  )

  # Build result
  result <- structure(
    merged,
    class = c("mutateR_consensus_batch", "list"),
    metadata = merged_meta
  )

  message("Merged results: ", merged_meta$n_successful, " total genes ",
          "(", primary_meta$n_successful, " primary + ",
          recovered_meta$n_successful, " recovered)")

  return(result)
}


#' Save recovery results
#'
#' Convenience function to save recovery outputs to an analysis directory.
#'
#' @param recovery_result Output from recover_failed_genes()
#' @param output_dir Analysis output directory
#' @param prefix File prefix (default "recovery")
#'
#' @return Paths to saved files (invisibly)
#'
#' @export
save_recovery_results <- function(recovery_result, output_dir, prefix = "recovery") {

  recovery_dir <- file.path(output_dir, "recovery")
  dir.create(recovery_dir, showWarnings = FALSE, recursive = TRUE)

  paths <- list()

  # Save recovered batch
  if (length(recovery_result$recovered) > 0) {
    paths$recovered <- file.path(recovery_dir, paste0(prefix, "_results.rds"))
    saveRDS(recovery_result$recovered, paths$recovered)
  }

  # Save still failed
  if (nrow(recovery_result$still_failed) > 0) {
    paths$still_failed <- file.path(recovery_dir, paste0(prefix, "_still_failed.rds"))
    saveRDS(recovery_result$still_failed, paths$still_failed)
  }

  # Save summary
  paths$summary <- file.path(recovery_dir, paste0(prefix, "_summary.rds"))
  saveRDS(recovery_result$summary, paths$summary)

  # Save excluded MHC genes if any
  if (length(recovery_result$excluded_mhc) > 0) {
    paths$mhc_excluded <- file.path(recovery_dir, paste0(prefix, "_mhc_excluded.rds"))
    saveRDS(recovery_result$excluded_mhc, paths$mhc_excluded)
  }

  message("Recovery results saved to: ", recovery_dir)
  invisible(paths)
}
