#' Venn diagram for proteomics analysis
#'
#' Generates a Venn diagram showing protein overlaps between groups.
#' Highlights proteins unique to each group and those shared across groups.
#' Based on presence/absence information from `proteo.normalize()`.
#'
#' @param normalized_data A list returned by `proteo.normalize()` (with remove_zeros = FALSE),
#'   containing the element `presence_per_group`, or a data frame with `presence_per_group`
#'   stored as an attribute.
#' @param cutoff Numeric. Minimum presence threshold to consider a protein present (default = 0).
#' @param results Logical. If TRUE, returns a list with proteins by group, unique proteins,
#'   and shared proteins (default = FALSE).
#' @param help Logical. If TRUE, prints this help message instead of running the function.
#'
#' @return Invisibly returns NULL, or if `results = TRUE`, a list with:
#'   \describe{
#'     \item{proteins_by_group}{List of proteins present in each group.}
#'     \item{proteins_unique}{List of unique proteins per group.}
#'     \item{proteins_shared}{Vector of proteins shared across all groups.}
#'   }
#'
#' @examples
#' \dontrun{
#' # --- Small example: presence/absence matrix ---
#' presence_per_group <- data.frame(
#'   ProteinID = paste0("P", 1:5),
#'   Control = c(1, 1, 0, 1, 0),
#'   Treatment = c(1, 0, 1, 1, 1)
#' )
#'
#' normalized_data <- list(presence_per_group = presence_per_group)
#'
#' # Generate Venn diagram and return results
#' venn_res <- proteo.venn(normalized_data, cutoff = 0, results = TRUE)
#' }
#'
#' @references
#' Venn, J. (1880). On the diagrammatic and mechanical representation of propositions and reasonings. \emph{The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science}, 10(59), 1–18. <doi:10.1080/14786448008627000>
#' Gentleman, R. et al. (2005). Bioconductor: open software development for computational biology and bioinformatics. \emph{Genome Biology}, 5(10), R80. <doi:10.1186/gb-2004-5-10-r80>
#'
#' @export
#' @importFrom ggvenn ggvenn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 theme_void theme element_blank

proteo.venn <- function(normalized_data, cutoff = 0, results = FALSE, help = FALSE) {

  # --- Help message ---
  if (help || missing(normalized_data)) {
    message("
Function proteo.venn()

Description:
  Generates a Venn diagram showing protein overlaps between groups.
  Highlights unique and shared proteins based on presence/absence data
  derived from proteo.normalize().

Usage:
  proteo.venn(normalized_data, cutoff = 0, results = FALSE)

Arguments:
  normalized_data  List returned by proteo.normalize() (presence_per_group) or data frame.
  cutoff           Minimum presence threshold to consider protein as present. Default 0.
  results          Logical. Return list with proteins by group, unique, and shared? Default FALSE.
  help             Logical. If TRUE, prints this help message.

Return:
  Invisibly returns NULL, or if results = TRUE, a list with:
    - proteins_by_group : List of proteins present in each group
    - proteins_unique   : List of unique proteins per group
    - proteins_shared   : Vector of proteins shared across all groups

Example:
  presence_per_group <- data.frame(
    ProteinID = paste0('P', 1:5),
    Control = c(1,1,0,1,0),
    Treatment = c(1,0,1,1,1)
  )
  normalized_data <- list(presence_per_group = presence_per_group)
  venn_res <- proteo.venn(normalized_data, cutoff=0, results=TRUE)
")
    return(invisible(NULL))
  }

  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    stop("Please install the 'ggvenn' package: install.packages('ggvenn')")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install the 'RColorBrewer' package: install.packages('RColorBrewer')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the 'ggplot2' package: install.packages('ggplot2')")
  }

  # --- Detect normalized object ---
  if (is.list(normalized_data) && !is.null(normalized_data$presence_per_group)) {
    presence_per_group <- normalized_data$presence_per_group
  } else if (!is.null(attr(normalized_data, "presence_per_group"))) {
    presence_per_group <- attr(normalized_data, "presence_per_group")
  } else {
    stop("The provided object does not contain 'presence_per_group'.
         Use proteo.normalize(..., remove_zeros = FALSE).")
  }

  # --- Ensure data.frame ---
  if (!is.data.frame(presence_per_group)) {
    presence_per_group <- as.data.frame(presence_per_group)
  }

  # --- Extract groups and proteins ---
  groups <- colnames(presence_per_group)[-1]  # exclude ProteinID
  protein_ids <- presence_per_group[[1]]      # first column = ProteinID

  # --- Proteins per group ---
  proteins_by_group <- lapply(groups, function(g) {
    protein_ids[presence_per_group[[g]] > cutoff]
  })
  names(proteins_by_group) <- groups

  # --- Unique proteins per group ---
  proteins_unique <- lapply(groups, function(g) {
    setdiff(proteins_by_group[[g]], unlist(proteins_by_group[groups != g]))
  })
  names(proteins_unique) <- groups

  # --- Shared proteins across all groups ---
  proteins_shared <- Reduce(intersect, proteins_by_group)

  # --- Venn diagram ---
  n_groups <- length(groups)
  n_colors <- max(3, n_groups)
  brewer_colors <- RColorBrewer::brewer.pal(n_colors, "Set1")[1:n_groups]

  p <- ggvenn::ggvenn(
    proteins_by_group,
    fill_color = brewer_colors,
    stroke_size = 0,
    set_name_size = 4,
    text_size = 4
  ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "right"
    )

  print(p)

  # --- Summary ---
  cat("\n=== Summary of proteins per group ===\n")
  for (grp in groups) {
    cat(sprintf("%s: %d proteins\n", grp, length(proteins_by_group[[grp]])))
  }
  cat("===============================\n\n")

  cat("=== Unique proteins per group ===\n")
  for (grp in groups) {
    cat(sprintf("%s: %d unique proteins\n", grp, length(proteins_unique[[grp]])))
  }
  cat(sprintf("Proteins shared across all groups: %d\n\n", length(proteins_shared)))

  if (!results) {
    message("Tip: To access detailed lists of proteins, call with 'results = TRUE'.\n")
  }

  if (results) {
    return(list(
      proteins_by_group = proteins_by_group,
      proteins_unique = proteins_unique,
      proteins_shared = proteins_shared
    ))
  } else {
    invisible(NULL)
  }
}










