#' Perform dimensionality reduction for proteomics data
#'
#' This function performs PCA and UMAP dimensionality reduction on proteomics data.
#' It handles normalized data from \code{proteo.normalize()}, ensures sample order matches
#' metadata, and replaces missing values with row means. Intra- and inter-group distances
#' are calculated on PC1 and PC2.
#'
#' @param normalized_data A data.frame or a list from \code{proteo.normalize()} containing
#'   normalized protein expression. Must have a ProteinID column.
#' @param metadata A data.frame containing sample metadata. Must contain columns for
#'   sample names and group assignment.
#' @param n_neighbors Integer. Number of neighbors for UMAP. If NULL or larger than the number
#'   of samples, defaults to half the number of samples. Default is NULL.
#' @param verbose Logical. If TRUE (default), prints messages and summary statistics.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{pca}{Result of \code{prcomp} PCA analysis.}
#'   \item{umap}{Result of \code{umap::umap} UMAP analysis.}
#'   \item{plot_pca}{ggplot2 PCA plot.}
#'   \item{plot_umap}{ggplot2 UMAP plot.}
#'   \item{n_neighbors}{Number of neighbors used in UMAP.}
#'   \item{dist_intra}{Mean intra-group distance (PC1 & PC2).}
#'   \item{dist_inter}{Mean inter-group distance (PC1 & PC2).}
#' }
#'
#' @examples
#' # Small example: raw data
#' raw_data <- data.frame(
#'   ProteinID = paste0("P", 1:5),
#'   Control1 = c(100, 200, 150, 300, 250),
#'   Control2 = c(110, 210, 160, 310, 260),
#'   Treatment1 = c(300, 100, 200, 150, 250),
#'   Treatment2 = c(310, 90, 210, 140, 260)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control1", "Control2", "Treatment1", "Treatment2"),
#'   Group = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' # Step 1: Normalize
#' normalized_obj <- proteo.normalize(raw_data, metadata)
#'
#' # Step 2: Perform dimensionality reduction
#' proteo.dimred(normalized_obj, metadata)
#'
#' @references
#' PCA: Pearson 1901 <doi:10.1080/14786440109462720>
#' UMAP: McInnes et al. 2020 <doi:10.21105/joss.02096>
#'
#' @export



proteo.dimred <- function(normalized_data, metadata,
                          n_neighbors = NULL, verbose = TRUE,
                          help = FALSE) {

  if (help || missing(normalized_data)) {
    message("
Function proteo.dimred()

Description:
  Performs PCA and UMAP dimensionality reduction on proteomics data.
  Handles normalized data from proteo.normalize(), fills missing values,
  and computes intra- and inter-group distances.

Usage:
  proteo.dimred(normalized_data, metadata, n_neighbors = NULL, verbose = TRUE)

Arguments:
  normalized_data  Proteomics object (proteoNorm) or data.frame with ProteinID.
  metadata         Data.frame with Sample and Group columns.
  n_neighbors      Number of neighbors for UMAP. Default NULL.
  verbose          Print messages and summary stats? Default TRUE.
  help             Show this message and example usage.

Example (functional):
  raw_data <- data.frame(
    ProteinID = paste0('P', 1:5),
    Control1   = c(100,200,150,300,250),
    Control2   = c(110,210,160,310,260),
    Treatment1 = c(300,100,200,150,250),
    Treatment2 = c(310,90,210,140,260)
  )

  metadata <- data.frame(
    Sample = c('Control1','Control2','Treatment1','Treatment2'),
    Group  = c('Control','Control','Treatment','Treatment')
  )

  # Step 1: Normalize data
  normalized_obj <- proteo.normalize(raw_data, metadata)

  # Step 2: Perform dimensionality reduction
  proteo.dimred(normalized_obj, metadata)
")
    return(invisible(NULL))
  }

  # --- Check required packages explicitly ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
  if (!requireNamespace("umap", quietly = TRUE)) stop("Please install 'umap'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")

  # --- Handle list input from proteo.normalize() ---
  if (is.list(normalized_data) && !is.data.frame(normalized_data)) {
    if ("log2_normalized" %in% names(normalized_data)) {
      normalized_data <- normalized_data[["log2_normalized"]]
    } else if ("log10_normalized" %in% names(normalized_data)) {
      normalized_data <- normalized_data[["log10_normalized"]]
    } else {
      stop("List object does not contain 'log2_normalized' or 'log10_normalized'.")
    }
  }

  # --- Detect ProteinID column ---
  possible_id_names <- c("ProteinID", "Protein_ID", "Protein", "Accession", "Protein.Accession", "protein_id")
  col_id <- which(names(normalized_data) %in% possible_id_names)
  if (length(col_id) != 1) stop("Unable to detect ProteinID column.")

  # --- Detect Sample and Group columns in metadata ---
  possible_sample_names <- c("Sample", "sample", "SAMPLE")
  col_sample <- which(names(metadata) %in% possible_sample_names)
  if (length(col_sample) != 1) stop("Unable to detect Sample column in metadata.")

  possible_group_names <- c("Group", "group", "GROUP", "Condition", "condition", "COND")
  col_group <- which(names(metadata) %in% possible_group_names)
  if (length(col_group) != 1) stop("Unable to detect Group column in metadata.")

  # --- Prepare expression matrix ---
  protein_names <- as.character(normalized_data[[col_id]])
  protein_names[is.na(protein_names) | protein_names == ""] <- "Protein_Unknown"
  protein_names <- make.unique(protein_names)

  expr_mat <- normalized_data[, -col_id, drop = FALSE]
  rownames(expr_mat) <- protein_names

  # --- Clean extreme/infinite values ---
  expr_mat[!is.finite(as.matrix(expr_mat))] <- NA
  expr_mat[expr_mat < -1e5 | expr_mat > 1e5] <- NA

  # --- Remove proteins with all NA ---
  expr_mat <- expr_mat[rowSums(is.na(expr_mat)) < ncol(expr_mat), , drop = FALSE]
  protein_names <- rownames(expr_mat)

  # --- Replace remaining NAs with row mean ---
  expr_mat <- t(apply(expr_mat, 1, function(x) {
    if (all(is.na(x))) return(rep(0, length(x)))
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  }))

  # --- Ensure sample order matches metadata ---
  samples <- metadata[[col_sample]]
  groups <- metadata[[col_group]]
  if (!all(colnames(expr_mat) == samples)) {
    expr_mat <- expr_mat[, samples, drop = FALSE]
  }

  n_samples <- ncol(expr_mat)

  # --- Adjust n_neighbors for UMAP ---
  if (is.null(n_neighbors) || n_neighbors >= n_samples) {
    n_neighbors <- max(2, floor(n_samples / 2))
    if(verbose) message(sprintf("n_neighbors for UMAP set to %d", n_neighbors))
  }

  # --- PCA ---
  pca_res <- stats::prcomp(t(expr_mat), center = TRUE, scale. = TRUE)
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

  df_pca <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = groups, Sample = samples)
  p_pca <- ggplot2::ggplot(df_pca, ggplot2::aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = "PCA - Dimensionality Reduction", x = "PC1", y = "PC2")

  # --- UMAP ---
  umap_res <- umap::umap(t(expr_mat), n_neighbors = n_neighbors, init = "random")
  df_umap <- data.frame(UMAP1 = umap_res$layout[,1], UMAP2 = umap_res$layout[,2], Group = groups, Sample = samples)
  p_umap <- ggplot2::ggplot(df_umap, ggplot2::aes(x = UMAP1, y = UMAP2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = "UMAP - Dimensionality Reduction", x = "UMAP1", y = "UMAP2")

  # --- Intra/inter-group distances (PC1 & PC2) ---
  dist_mat <- as.matrix(stats::dist(df_pca[,c("PC1","PC2")]))
  intra_dists <- c(); inter_dists <- c()
  for(i in 1:(n_samples-1)) {
    for(j in (i+1):n_samples) {
      if(df_pca$Group[i] == df_pca$Group[j]) {
        intra_dists <- c(intra_dists, dist_mat[i,j])
      } else {
        inter_dists <- c(inter_dists, dist_mat[i,j])
      }
    }
  }

  mean_intra <- mean(intra_dists)
  mean_inter <- mean(inter_dists)

  # --- Report ---
  if (verbose) {
    report <- paste0(
      "\n=== Dimensionality Reduction Report ===\n",
      sprintf("Variance explained PC1: %.1f%%\n", var_explained[1] * 100),
      sprintf("Variance explained PC2: %.1f%%\n", var_explained[2] * 100),
      sprintf("UMAP n_neighbors: %d\n", n_neighbors),
      sprintf("Mean intra-group distance (PC1 & PC2): %.3f\n", mean_intra),
      sprintf("Mean inter-group distance (PC1 & PC2): %.3f\n", mean_inter),
      "=========================\n"
    )
    message(report)
  }

  # --- Plot ---
  print(p_pca)
  print(p_umap)

  invisible(list(
    pca = pca_res,
    umap = umap_res,
    plot_pca = p_pca,
    plot_umap = p_umap,
    n_neighbors = n_neighbors,
    dist_intra = mean_intra,
    dist_inter = mean_inter
  ))
}





