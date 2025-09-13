#' Volcano Heatmap of Proteomics Data
#'
#' Creates a heatmap of the top most variable proteins across samples or groups.
#' Automatically applies per-protein z-scoring, selects the most variable proteins,
#' and annotates samples by group. Returns the selected proteins and data if requested.
#'
#' @param normalized_data A `data.frame` or `proteoNorm` object containing normalized protein abundances.
#'   Must include a column named `ProteinID`.
#' @param metadata A `data.frame` with at least two columns:
#'   - `Sample`: matching the column names in `normalized_data`
#'   - `Group`: group assignment for each sample
#' @param group_colname Column name in `metadata` containing group information. Default = `"Group"`.
#' @param top_n Integer. Number of most variable proteins to include. Default = 100.
#' @param cluster_rows Logical. Whether to cluster rows in the heatmap. Default = TRUE.
#' @param cluster_cols Logical. Whether to cluster columns in the heatmap. Default = TRUE.
#' @param show_rownames Logical. Show protein IDs in the heatmap. Default = FALSE.
#' @param show_colnames Logical. Show sample names in the heatmap. Default = TRUE.
#' @param palette Color palette for the heatmap. Default = blue-white-red.
#' @param results Logical. If TRUE, returns a list with selected proteins and data. Default = FALSE.
#' @param assign_result Logical. If TRUE, assigns the result to the environment. Default = FALSE.
#' @param assign_name Name of the object to assign if `assign_result = TRUE`. Default = `"heatmap_data"`.
#' @param envir Environment where the result will be assigned. Default = `parent.frame()`.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return
#' If `results = TRUE`, returns a list with:
#' \describe{
#'   \item{top_proteins}{Vector of top selected proteins.}
#'   \item{expr_z_score}{Z-scored expression matrix.}
#'   \item{annotation_col}{Sample annotations for the heatmap.}
#' }
#' If `assign_result = TRUE`, the same list is assigned to the chosen environment.
#' Otherwise, only the heatmap is displayed.
#'
#' @examples
#' \dontrun{
#' # Small example dataset
#' raw_data <- data.frame(
#'   ProteinID = c('P1','P2','P3','P4'),
#'   Control_1 = c(6.2,0,4.5,3.1),
#'   Control_2 = c(6.1,0,4.8,3.3),
#'   Treatment_1 = c(5.9,0,5.0,2.9),
#'   Treatment_2 = c(6.0,0,4.9,3.0)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c('Control_1','Control_2','Treatment_1','Treatment_2'),
#'   Group  = c('Control','Control','Treatment','Treatment')
#' )
#'
#' # Step 1: Normalize data
#' normalized_obj <- proteo.normalize(raw_data, metadata)
#'
#' # Step 2: Generate heatmap
#' proteo.heatmap(normalized_obj, metadata, top_n = 3)
#' proteo.heatmap(normalized_obj, metadata, top_n = 3, results = TRUE)
#' }
#'
#' @references
#' Kolde R (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
#'
#' @export

proteo.heatmap <- function(normalized_data, metadata,
                           group_colname = "Group",
                           top_n = 100,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           show_rownames = FALSE,
                           show_colnames = TRUE,
                           palette = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
                           results = FALSE,
                           assign_result = FALSE,
                           assign_name = "heatmap_data",
                           envir = parent.frame(),
                           help = FALSE) {

  if (help || missing(normalized_data) || missing(metadata)) {
    message("
Function proteo.heatmap()

Description:
  Creates a heatmap of the top most variable proteins across samples or groups.
  Applies z-score normalization per protein and annotates samples by group.
  Returns heatmap and, optionally, the processed data.

Usage:
  proteo.heatmap(normalized_data, metadata, top_n = 100, ...)

Arguments:
  normalized_data  Normalized data.frame or proteoNorm object with ProteinID column.
  metadata         Data frame with 'Sample' and group assignment columns.
  top_n            Number of most variable proteins to include (default 100).
  cluster_rows     Cluster rows? Default TRUE.
  cluster_cols     Cluster columns? Default TRUE.
  show_rownames    Show protein IDs? Default FALSE.
  show_colnames    Show sample names? Default TRUE.
  palette          Color palette. Default blue-white-red.
  results          Return processed data? Default FALSE.
  assign_result    Assign result to environment? Default FALSE.
  assign_name      Name of assigned object. Default 'heatmap_data'.
  envir            Environment to assign if assign_result = TRUE.
  help             Show this message.

Example (functional):
  raw_data <- data.frame(
    ProteinID = c('P1','P2','P3','P4'),
    Control_1 = c(6.2,0,4.5,3.1),
    Control_2 = c(6.1,0,4.8,3.3),
    Treatment_1 = c(5.9,0,5.0,2.9),
    Treatment_2 = c(6.0,0,4.9,3.0)
  )

  metadata <- data.frame(
    Sample = c('Control_1','Control_2','Treatment_1','Treatment_2'),
    Group  = c('Control','Control','Treatment','Treatment')
  )

  # Step 1: Normalize data
  normalized_obj <- proteo.normalize(raw_data, metadata)

  # Step 2: Generate heatmap
  proteo.heatmap(normalized_obj, metadata, top_n = 3)
")
    return(invisible(NULL))
  }

  # --- Required packages ---
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install 'pheatmap'.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")

  # --- Extract normalized dataframe ---
  if ("proteoNorm" %in% class(normalized_data)) {
    df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
    if (length(df_names) != 1) stop("Could not detect the normalized data.frame in 'proteoNorm'.")
    normalized_data <- normalized_data[[df_names]]
  }

  # --- Prepare numeric matrix ---
  if (!"ProteinID" %in% colnames(normalized_data)) stop("Missing column 'ProteinID' in normalized_data.")
  protein_ids <- normalized_data$ProteinID
  expr_data <- normalized_data[, -which(colnames(normalized_data) == "ProteinID")]
  expr_data <- as.matrix(expr_data)
  storage.mode(expr_data) <- "numeric"
  rownames(expr_data) <- protein_ids

  # --- Prepare groups ---
  if (!"Sample" %in% colnames(metadata)) stop("Metadata must contain the column 'Sample'.")
  if (!(group_colname %in% colnames(metadata))) stop("Group column not found in metadata.")

  sample_groups <- metadata[[group_colname]]
  names(sample_groups) <- metadata$Sample

  # --- Z-score per protein ---
  expr_z <- t(scale(t(expr_data)))
  message("Applied z-score per protein.")

  # --- Select top proteins by variance ---
  var_proteins <- apply(expr_z, 1, stats::var, na.rm = TRUE)
  top_proteins <- names(sort(var_proteins, decreasing = TRUE))[1:min(top_n, length(var_proteins))]
  expr_z_top <- expr_z[top_proteins, , drop = FALSE]
  message(sprintf("Selected %d proteins with highest variance.", length(top_proteins)))

  # --- Sample annotation ---
  annotation_col <- data.frame(Group = sample_groups[colnames(expr_z_top)])
  rownames(annotation_col) <- colnames(expr_z_top)

  # --- Generate heatmap ---
  pheatmap::pheatmap(expr_z_top,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     annotation_col = annotation_col,
                     color = palette,
                     border_color = NA,
                     main = sprintf("Top %d Most Variable Proteins", length(top_proteins)),
                     fontsize = 10)

  # --- Result object ---
  result <- list(
    top_proteins = top_proteins,
    expr_z_score = expr_z_top,
    annotation_col = annotation_col
  )

  if (assign_result) {
    assign(assign_name, result, envir = envir)
    message(sprintf("Object '%s' created in the current environment.", assign_name))
  }

  if (results) {
    return(result)
  } else {
    message("Tip: To inspect the selected proteins and z-score matrix, call with 'results = TRUE'")
    invisible(NULL)
  }
}




