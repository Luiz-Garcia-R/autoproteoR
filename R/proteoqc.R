#' Quality Control Plots for Normalized Proteomics Data
#'
#' Generates quality control (QC) plots for normalized proteomics data.
#' Includes a boxplot of intensity distributions per sample and a PCA plot colored by group.
#' The input must be an object returned by \code{proteo.normalize()} (class \code{proteoNorm}).
#'
#' @param normalized_data Object of class \code{proteoNorm} returned by \code{proteo.normalize()}.
#' @param metadata A data frame containing sample metadata. Must include:
#'                 - \code{Sample}: column matching the column names of the normalized data (excluding ProteinID)
#'                 - \code{Group}: group/condition labels for PCA coloring
#' @param cluster_rows Logical. Whether to cluster rows in the boxplot? Default FALSE.
#' @param cluster_cols Logical. Whether to cluster columns in the boxplot? Default FALSE.
#' @param angle_col Numeric. Angle of column labels in the boxplot. Default 45.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return A list with:
#'   \item{boxplot}{A ggplot2 object showing distribution of intensities per sample}
#'   \item{pca}{A ggplot2 object showing PCA colored by group}
#' @examples
#' # Example workflow
#' raw_data <- data.frame(
#'   ProteinID = c("P1","P2","P3","P4"),
#'   Control_1 = c(6.2,0,4.5,3.1),
#'   Control_2 = c(6.1,0,4.8,3.3),
#'   Treatment_1 = c(5.9,0,5.0,2.9),
#'   Treatment_2 = c(6.0,0,4.9,3.0)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control_1","Control_2","Treatment_1","Treatment_2"),
#'   Group  = c("Control","Control","Treatment","Treatment")
#' )
#'
#' # Step 1: Normalize
#' normalized_obj <- proteo.normalize(raw_data, metadata)
#'
#' # Step 2: QC plots
#' proteo.qc(normalized_obj, metadata)
#'
#' @references
#' Bolstad, B. M., Irizarry, R. A., Åstrand, M., & Speed, T. P. (2003). A comparison of normalization methods for high density oligonucleotide array data. \emph{Bioinformatics}, 19(2), 185–193. <doi:10.1093/bioinformatics/btg405>
#' Jolliffe, I. T. (2002). Principal Component Analysis. Springer Series in Statistics. <doi:10.1007/b98835>
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @export


proteo.qc <- function(normalized_data, metadata,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      angle_col = 45, help = FALSE) {

  # --- Help message ---
  if (help || missing(normalized_data) || missing(metadata)) {
    message("
Function proteo.qc()

Description:
  Generates quality control (QC) plots for normalized proteomics data.
  Includes a boxplot of intensity distributions per sample and a PCA plot
  colored by group. Helps assess data quality, detect potential outliers,
  and visualize data structure.

Important:
  The input 'normalized_data' should be an object returned by proteo.normalize()
  (class 'proteoNorm') containing a <method>_normalized data.frame.
  Directly using a raw data.frame will produce an error.

Usage:
  proteo.qc(normalized_data, metadata, cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45)

Arguments:
  normalized_data  Object of class 'proteoNorm' (from proteo.normalize()).
  metadata         Data frame with sample metadata. Must include:
                     - 'Sample': matches column names of normalized_data (excluding 'ProteinID')
                     - 'Group' : group/condition labels for PCA coloring
  cluster_rows     Logical. Whether to cluster rows in boxplot? Default FALSE.
  cluster_cols     Logical. Whether to cluster columns in boxplot? Default FALSE.
  angle_col        Numeric. Angle of column labels in boxplot. Default 45.
  help             Logical. If TRUE, prints this help message.

Return:
  A list with:
    - boxplot  : ggplot2 object showing distribution of intensities per sample
    - pca      : ggplot2 object showing PCA colored by group

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

  # Step 2: Generate QC plots
  proteo.qc(normalized_obj, metadata)
")
    return(invisible(NULL))
  }

  # ---- check packages ----
  req_pkgs <- c("ggplot2", "pheatmap", "tidyr")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  # ---- extract normalized data.frame from proteoNorm/list if needed ----
  expr_df <- NULL
  if (is.list(normalized_data)) {
    # prefer element matching *_normalized
    df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
    if (length(df_names) >= 1) {
      expr_df <- normalized_data[[df_names[1]]]
    } else {
      # fallback: first element that is a data.frame and contains ProteinID
      df_candidates <- Filter(function(x) is.data.frame(x) && "ProteinID" %in% colnames(x), normalized_data)
      if (length(df_candidates) >= 1) expr_df <- df_candidates[[1]]
    }
  } else if (is.data.frame(normalized_data)) {
    expr_df <- normalized_data
  }

  if (is.null(expr_df) || !is.data.frame(expr_df)) {
    stop("Input must be a data.frame or a list returned by proteo.normalize() (containing a <name>_normalized data.frame).")
  }

  # ensure not a tibble with rownames weirdness
  expr_df <- as.data.frame(expr_df, stringsAsFactors = FALSE, check.names = FALSE)

  # ---- checks ----
  if (!"ProteinID" %in% colnames(expr_df)) stop("The normalized data.frame must contain a column named 'ProteinID'.")
  if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame with a column named 'Sample'.")

  if (!"Sample" %in% colnames(metadata)) stop("`metadata` must contain a column named 'Sample'.")

  # ---- build expression matrix, align samples ----
  data_mat <- expr_df[, setdiff(colnames(expr_df), "ProteinID"), drop = FALSE]
  # keep only numeric columns (safety)
  numeric_cols <- vapply(data_mat, is.numeric, logical(1))
  if (!all(numeric_cols)) {
    warning("Some non-numeric columns found and will be dropped from expression matrix: ",
            paste(names(numeric_cols)[!numeric_cols], collapse = ", "))
    data_mat <- data_mat[, numeric_cols, drop = FALSE]
  }
  rownames(data_mat) <- as.character(expr_df$ProteinID)

  # match metadata samples to expression columns
  common_samples <- intersect(metadata$Sample, colnames(data_mat))
  if (length(common_samples) == 0) stop("No sample names match between metadata$Sample and expression matrix column names.")
  # reorder columns to metadata order (only keep matching)
  meta_idx <- match(common_samples, metadata$Sample)
  metadata_sub <- metadata[meta_idx, , drop = FALSE]
  data_mat <- as.data.frame(data_mat[, common_samples, drop = FALSE], check.names = FALSE)

  # ---- correlation matrix across all samples ----
  corr_mat_all <- stats::cor(as.matrix(data_mat), use = "pairwise.complete.obs", method = "pearson")

  # ---- matrix of "r = ..." strings for display ----
  corr_text <- matrix(paste0("r = ", formatC(round(corr_mat_all, 2), format = "f", digits = 2)),
                      nrow = nrow(corr_mat_all),
                      ncol = ncol(corr_mat_all),
                      dimnames = dimnames(corr_mat_all))

  # ---- long dataframe for ggplots ----
  df_long <- data.frame(ProteinID = rownames(data_mat), data_mat, check.names = FALSE, stringsAsFactors = FALSE)
  df_long <- tidyr::pivot_longer(df_long, cols = -ProteinID, names_to = "Sample", values_to = "Expression")
  # attach metadata (left join)
  df_long <- merge(df_long, metadata_sub, by = "Sample", all.x = TRUE, sort = FALSE)

  # ---- density plot ----
  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(x = Expression, color = Sample, fill = Sample)) +
    ggplot2::geom_density(alpha = 0.3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Expression distribution per sample",
                  x = "Normalized expression", y = "Density") +
    ggplot2::theme(legend.position = "bottom")

  # ---- boxplot ----
  p_boxplot <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = Expression, fill = Sample)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Expression boxplot per sample", x = "", y = "Normalized expression") +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = angle_col, hjust = 1))

  print(p_density)
  print(p_boxplot)

  # ---- correlation heatmap ----
  pheatmap::pheatmap(
    corr_mat_all,
    main = "Correlation heatmap across all samples",
    silent = FALSE,
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    clustering_method = "complete",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    fontsize = 10,
    legend = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = ifelse(cluster_rows, 30, 0),
    treeheight_col = ifelse(cluster_cols, 30, 0),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    display_numbers = corr_text,
    number_color = "black",
    angle_col = angle_col
  )

  # ---- QC summary ----
  n_proteins <- nrow(data_mat)
  n_samples <- ncol(data_mat)
  mean_global <- suppressWarnings(mean(as.matrix(data_mat), na.rm = TRUE))
  mean_corr_all <- mean(corr_mat_all[lower.tri(corr_mat_all)], na.rm = TRUE)

  message("=== QC Summary ===")
  cat(sprintf("Number of proteins analyzed: %d\n", n_proteins))
  cat(sprintf("Number of samples: %d\n", n_samples))
  cat(sprintf("Global mean of normalized expression: %.3f\n", mean_global))
  cat(sprintf("Average correlation between all samples: %.3f\n", mean_corr_all))
  if (mean_corr_all < 0.7) {
    message("Warning: Low average correlation may indicate variability or normalization issues.")
  }
  message("==================")

  # ---- return ----
  invisible(list(
    data_ready = data_mat,
    protein_ids = rownames(data_mat),
    corr_matrix_all = corr_mat_all,
    plot_density = p_density,
    plot_boxplot = p_boxplot
  ))
}


