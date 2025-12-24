#' Score gRNAs using multiple on-target scoring methods
#'
#' Runs multiple scoring methods on the same set of gRNAs and returns
#' all scores in a single data frame for comparative analysis.
#'
#' @param grna_gr GRanges returned by find_cas9_sites() or find_cas12a_sites().
#' @param methods Character vector of scoring methods to run.
#'        For Cas9: any subset of c("ruleset1", "deepspcas9", "ruleset3", "deephf")
#'        For Cas12a: any subset of c("deepcpf1", "enpamgb")
#'        Default NULL auto-selects based on detected nuclease.
#' @param tracr Character. tracrRNA for RuleSet3 (default "Chen2013").
#' @param deephf_var Character. DeepHF variant (default "wt_u6").
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return A data.frame with columns:
#'   - grna_id: Unique identifier for each gRNA
#'   - protospacer_sequence: The protospacer sequence
#'   - exon_rank: Exon location
#'   - One column per scoring method containing raw scores
#'
#' @examples
#' \dontrun
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Get gRNA sites
#' exons <- get_exon_structures("ENST00000269305", "hsapiens", output = "GRanges")
#' sites <- find_cas9_sites(exons, BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Score with multiple methods
#' multi_scores <- score_grnas_multi(sites,
#'                                   methods = c("ruleset1", "deepspcas9", "ruleset3", "deephf"))
#'
#' # Explore correlations
#' plot_score_correlations(multi_scores)
#' }
#'
#' @export
score_grnas_multi <- function(grna_gr,
                              methods = NULL,
                              tracr = "Chen2013",
                              deephf_var = "wt_u6",
                              quiet = FALSE) {

  if (!inherits(grna_gr, "GRanges")) {
    stop("Input must be a GRanges object from find_cas9_sites() or find_cas12a_sites().")
  }

  if (length(grna_gr) == 0) {
    warning("Empty GRanges provided.")
    return(data.frame())
  }

  # --- Detect nuclease type from sequence context length ---
  ctx_len <- nchar(mcols(grna_gr)$sequence_context[1])

  if (ctx_len == 30) {
    nuclease <- "Cas9"
    available_methods <- c("ruleset1", "deepspcas9", "ruleset3", "deephf")
  } else if (ctx_len == 34) {
    nuclease <- "Cas12a"
    available_methods <- c("deepcpf1", "enpamgb")
  } else {
    stop("Unrecognised sequence context length: ", ctx_len,
         ". Expected 30 (Cas9) or 34 (Cas12a).")
  }

  # --- Set default methods if not specified ---
  if (is.null(methods)) {
    methods <- available_methods
    if (!quiet) message("Auto-selected methods for ", nuclease, ": ",
                        paste(methods, collapse = ", "))
  }

  # --- Validate requested methods ---
  invalid <- setdiff(methods, available_methods)
  if (length(invalid) > 0) {
    stop("Invalid methods for ", nuclease, ": ", paste(invalid, collapse = ", "),
         "\nAvailable: ", paste(available_methods, collapse = ", "))
  }

  # --- Build base data frame ---
  n_grnas <- length(grna_gr)

  result_df <- data.frame(
    grna_id = paste0("gRNA_", seq_len(n_grnas)),
    protospacer_sequence = as.character(mcols(grna_gr)$protospacer_sequence),
    sequence_context = as.character(mcols(grna_gr)$sequence_context),
    exon_rank = mcols(grna_gr)$exon_rank,
    stringsAsFactors = FALSE
  )

  # --- Add genomic coordinates ---
  result_df$chr <- as.character(seqnames(grna_gr))
  result_df$start <- start(grna_gr)
  result_df$end <- end(grna_gr)
  result_df$strand <- as.character(strand(grna_gr))

  # --- Run each scoring method ---
  for (method in methods) {
    if (!quiet) message("Scoring with ", method, "...")

    scored_gr <- tryCatch({
      score_grnas(grna_gr,
                  method = method,
                  tracr = tracr,
                  deephf_var = deephf_var)
    }, error = function(e) {
      warning("Method '", method, "' failed: ", e$message)
      return(NULL)
    })

    if (!is.null(scored_gr)) {
      # Extract scores and add to result
      scores <- as.numeric(mcols(scored_gr)$ontarget_score)
      result_df[[method]] <- scores
    } else {
      result_df[[method]] <- NA_real_
    }
  }

  # --- Add metadata attributes ---
  attr(result_df, "nuclease") <- nuclease
  attr(result_df, "methods") <- methods
  attr(result_df, "n_grnas") <- n_grnas

  # --- Add method classification ---
  if (nuclease == "Cas9") {
    attr(result_df, "modern_methods") <- intersect(c("deepspcas9", "deephf", "ruleset3"), methods)
    attr(result_df, "legacy_methods") <- intersect(c("ruleset1"), methods)
  } else {
    attr(result_df, "modern_methods") <- methods
    attr(result_df, "legacy_methods") <- character(0)
  }

  if (!quiet) {
    message("Completed scoring ", n_grnas, " gRNAs with ",
            length(methods), " methods.")
  }

  return(result_df)
}
