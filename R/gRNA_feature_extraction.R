#' @title gRNA Feature Extraction
#' @description
#' Computes sequence-based features for gRNA SHAP analysis.
#' Features include compositional, positional, structural, and thermodynamic properties.
#'
#' @name gRNA_feature_extraction
NULL


#' Extract sequence features for SHAP analysis
#'
#' Computes a comprehensive set of sequence-based features for each gRNA
#' to enable Shapley value analysis of scoring method behavior.
#'
#' @param scores_df Data frame from score_grnas_multi()
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{grna_id}{gRNA identifier from input}
#'     \item{gc_content}{Overall GC fraction (0-1)}
#'     \item{gc_count}{Absolute GC count}
#'     \item{gc_proximal}{GC fraction in PAM-proximal 10nt (positions 11-20)}
#'     \item{gc_distal}{GC fraction in PAM-distal 10nt (positions 1-10)}
#'     \item{max_homopolymer}{Length of longest homopolymer run}
#'     \item{has_polyT/A/G/C}{Binary flags for NNNN runs}
#'     \item{freq_XX}{Dinucleotide frequencies (16 features)}
#'     \item{posN_X}{Position-specific nucleotide indicators (80 features)}
#'     \item{estimated_tm}{Wallace rule Tm approximation}
#'     \item{self_complement_score}{Simple self-complementarity metric}
#'     \item{linguistic_complexity}{K-mer diversity measure}
#'     \item{pre_pam_X}{Nucleotide at position 20 (4 binary features)}
#'     \item{exon_rank}{Exon position from input}
#'   }
#'
#' @details
#' Features are designed for Cas9 with 20nt protospacers. The function
#' automatically handles sequences of different lengths but interpretation
#' of positional features assumes standard Cas9 geometry.
#'
#' Feature categories:
#' \describe{
#'   \item{Compositional}{GC content (global, regional), nucleotide counts}
#'   \item{Positional}{Position-specific nucleotide indicators for all 20 positions}
#'   \item{Sequence patterns}{Dinucleotide frequencies, homopolymer runs}
#'   \item{Structural proxies}{Self-complementarity, linguistic complexity}
#'   \item{Thermodynamic proxies}{Estimated Tm (Wallace rule)}
#' }
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Get scored gRNAs
#' exons <- get_exon_structures("ENST00000269305", "hsapiens", output = "GRanges")
#' sites <- find_cas9_sites(exons, BSgenome.Hsapiens.UCSC.hg38)
#' scores <- score_grnas_multi(sites)
#'
#' # Extract features
#' features <- extract_grna_features(scores)
#' dim(features)  # n_gRNAs x 113 features
#' }
#'
#' @seealso \code{\link{compute_shap_targets}} for target variable computation,
#'   \code{\link{analyze_feature_importance}} for SHAP analysis
#'
#' @export
extract_grna_features <- function(scores_df) {


  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' required. Install with: install.packages('stringr')")
  }

  # Validate input

  if (!"protospacer_sequence" %in% names(scores_df)) {
    stop("scores_df must contain 'protospacer_sequence' column")
  }

  seqs <- toupper(scores_df$protospacer_sequence)
  n <- length(seqs)

  # Initialize output
  features <- data.frame(
    grna_id = if ("grna_id" %in% names(scores_df)) scores_df$grna_id else paste0("gRNA_", seq_len(n)),
    stringsAsFactors = FALSE
  )

  # ═══════════════════════════════════════════════════════════════
  # 1. Global GC Composition
  # ═══════════════════════════════════════════════════════════════

  features$gc_content <- vapply(seqs, function(s) {
    (stringr::str_count(s, "G") + stringr::str_count(s, "C")) / nchar(s)
  }, numeric(1))

  features$gc_count <- vapply(seqs, function(s) {
    stringr::str_count(s, "G") + stringr::str_count(s, "C")
  }, numeric(1))

  # ═══════════════════════════════════════════════════════════════
  # 2. Regional GC (PAM-proximal vs PAM-distal)
  # ═══════════════════════════════════════════════════════════════

  features$gc_proximal <- vapply(seqs, function(s) {
    # Last 10 nt (PAM-proximal for Cas9, positions 11-20)
    seq_len <- nchar(s)
    if (seq_len < 10) return(NA_real_)
    prox <- substr(s, seq_len - 9, seq_len)
    (stringr::str_count(prox, "G") + stringr::str_count(prox, "C")) / 10
  }, numeric(1))

  features$gc_distal <- vapply(seqs, function(s) {
    # First 10 nt (PAM-distal, positions 1-10)
    if (nchar(s) < 10) return(NA_real_)
    dist <- substr(s, 1, 10)
    (stringr::str_count(dist, "G") + stringr::str_count(dist, "C")) / 10
  }, numeric(1))

  # ═══════════════════════════════════════════════════════════════
  # 3. Homopolymer Runs
  # ═══════════════════════════════════════════════════════════════

  features$max_homopolymer <- vapply(seqs, function(s) {
    runs <- gregexpr("(.)\\1+", s, perl = TRUE)[[1]]
    if (runs[1] == -1) return(1L)
    max(attr(runs, "match.length"))
  }, integer(1))

  features$has_polyT <- as.integer(grepl("TTTT", seqs))
  features$has_polyA <- as.integer(grepl("AAAA", seqs))
  features$has_polyG <- as.integer(grepl("GGGG", seqs))
  features$has_polyC <- as.integer(grepl("CCCC", seqs))

  # ═══════════════════════════════════════════════════════════════
  # 4. Dinucleotide Frequencies
  # ═══════════════════════════════════════════════════════════════

  dinucs <- c("AA", "AT", "AG", "AC",
              "TA", "TT", "TG", "TC",
              "GA", "GT", "GG", "GC",
              "CA", "CT", "CG", "CC")

  for (di in dinucs) {
    features[[paste0("freq_", di)]] <- vapply(seqs, function(s) {
      stringr::str_count(s, di) / (nchar(s) - 1)
    }, numeric(1))
  }

  # ═══════════════════════════════════════════════════════════════
  # 5. Position-Specific Nucleotides
  # ═══════════════════════════════════════════════════════════════

  # All 20 positions for Cas9
  for (pos in 1:20) {
    for (nuc in c("A", "T", "G", "C")) {
      features[[paste0("pos", pos, "_", nuc)]] <- vapply(seqs, function(s) {
        if (nchar(s) < pos) return(NA_integer_)
        as.integer(substr(s, pos, pos) == nuc)
      }, integer(1))
    }
  }

  # ═══════════════════════════════════════════════════════════════
  # 6. Thermodynamic Proxies
  # ═══════════════════════════════════════════════════════════════

  # Wallace rule approximation for Tm
  # Tm = 2(A+T) + 4(G+C) for oligos < 14 nt
  # For longer oligos this is a rough proxy
  features$estimated_tm <- vapply(seqs, function(s) {
    gc <- stringr::str_count(s, "[GC]")
    at <- stringr::str_count(s, "[AT]")
    2 * at + 4 * gc
  }, numeric(1))

  # ═══════════════════════════════════════════════════════════════
  # 7. Self-Complementarity Score
  # ═══════════════════════════════════════════════════════════════

  features$self_complement_score <- vapply(seqs, function(s) {
    # Count matching dinucleotides between sequence and its reverse complement
    rc <- chartr("ATGC", "TACG", s)
    rc <- paste(rev(strsplit(rc, "")[[1]]), collapse = "")
    sum(vapply(1:(nchar(s) - 1), function(i) {
      as.integer(substr(s, i, i + 1) == substr(rc, i, i + 1))
    }, integer(1)))
  }, numeric(1))

  # ═══════════════════════════════════════════════════════════════
  # 8. Linguistic Complexity
  # ═══════════════════════════════════════════════════════════════

  features$linguistic_complexity <- vapply(seqs, function(s) {
    n_chars <- nchar(s)
    if (n_chars < 3) return(NA_real_)

    # Count unique k-mers vs theoretical maximum
    kmers_2 <- vapply(1:(n_chars - 1), function(i) substr(s, i, i + 1), character(1))
    kmers_3 <- vapply(1:(n_chars - 2), function(i) substr(s, i, i + 2), character(1))

    observed <- length(unique(kmers_2)) + length(unique(kmers_3))
    expected <- min(16, n_chars - 1) + min(64, n_chars - 2)

    observed / expected
  }, numeric(1))

  # ═══════════════════════════════════════════════════════════════
  # 9. PAM-Adjacent Context (Position 20 = pre-PAM)
  # ═══════════════════════════════════════════════════════════════

  pre_pam_nuc <- vapply(seqs, function(s) {
    substr(s, nchar(s), nchar(s))
  }, character(1))

  for (nuc in c("A", "T", "G", "C")) {
    features[[paste0("pre_pam_", nuc)]] <- as.integer(pre_pam_nuc == nuc)
  }

  # ═══════════════════════════════════════════════════════════════
  # 10. Exon Position (from input)
  # ═══════════════════════════════════════════════════════════════

  if ("exon_rank" %in% names(scores_df)) {
    features$exon_rank <- scores_df$exon_rank
  }

  return(features)
}


#' Get list of all available feature names
#'
#' Returns the names of all features computed by \code{extract_grna_features()}.
#' Useful for feature selection and ablation studies.
#'
#' @param include_id Logical. Include 'grna_id' column (default FALSE)
#'
#' @return Character vector of feature names
#'
#' @examples
#' feature_names <- get_feature_names()
#' length(feature_names)
#'
#' # Categorize features
#' positional <- grep("^pos\\d+_", feature_names, value = TRUE)
#' dinuc <- grep("^freq_", feature_names, value = TRUE)
#'
#' @export
get_feature_names <- function(include_id = FALSE) {

  # Build feature list programmatically
  features <- c(
    # GC composition
    "gc_content", "gc_count", "gc_proximal", "gc_distal",

    # Homopolymers
    "max_homopolymer", "has_polyT", "has_polyA", "has_polyG", "has_polyC"
  )

  # Dinucleotides
  dinucs <- c("AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC",
              "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC")
  features <- c(features, paste0("freq_", dinucs))

  # Positional (20 positions x 4 nucleotides)
  for (pos in 1:20) {
    for (nuc in c("A", "T", "G", "C")) {
      features <- c(features, paste0("pos", pos, "_", nuc))
    }
  }

  # Thermodynamic and structural
  features <- c(features,
                "estimated_tm",
                "self_complement_score",
                "linguistic_complexity")

  # PAM-adjacent
  features <- c(features,
                "pre_pam_A", "pre_pam_T", "pre_pam_G", "pre_pam_C")

  # Exon position
  features <- c(features, "exon_rank")

  if (include_id) {
    features <- c("grna_id", features)
  }

  return(features)
}


#' Categorize features by type
#'
#' Returns a named list grouping feature names by their category.
#' Useful for understanding feature importance results.
#'
#' @return Named list with categories:
#'   \describe{
#'     \item{composition}{GC content features}
#'     \item{homopolymer}{Homopolymer run features}
#'     \item{dinucleotide}{Dinucleotide frequency features}
#'     \item{positional}{Position-specific nucleotide features}
#'     \item{thermodynamic}{Tm and energy-related features}
#'     \item{structural}{Self-complementarity, complexity features}
#'     \item{context}{PAM-adjacent and exon features}
#'   }
#'
#' @examples
#' cats <- get_feature_categories()
#' names(cats)
#' length(cats$positional)
#'
#' @export
get_feature_categories <- function() {

  dinucs <- c("AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC",
              "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC")

  positional <- character()
  for (pos in 1:20) {
    for (nuc in c("A", "T", "G", "C")) {
      positional <- c(positional, paste0("pos", pos, "_", nuc))
    }
  }

  list(
    composition = c("gc_content", "gc_count", "gc_proximal", "gc_distal"),
    homopolymer = c("max_homopolymer", "has_polyT", "has_polyA", "has_polyG", "has_polyC"),
    dinucleotide = paste0("freq_", dinucs),
    positional = positional,
    thermodynamic = c("estimated_tm"),
    structural = c("self_complement_score", "linguistic_complexity"),
    context = c("pre_pam_A", "pre_pam_T", "pre_pam_G", "pre_pam_C", "exon_rank")
  )
}
