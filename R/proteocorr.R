#' Compute and visualize correlation of proteomics data
#'
#' Computes pairwise correlations among samples or groups for proteomics data.
#' Automatically chooses Pearson or Spearman correlation based on normality and ties,
#' and provides heatmap and scatterplot visualizations.
#'
#' @param normalized_data A \code{data.frame} or \code{proteoNorm} object containing normalized protein expression.
#' @param metadata Optional \code{data.frame} with sample metadata. Required if \code{group_level = TRUE}.
#' @param group_level Logical. If TRUE (default), correlation is computed at the group level (averaging samples per group).
#' @param method Character. Correlation method: "pearson", "spearman", "kendall", or "auto" (default).
#' @param plot Logical. If TRUE (default), a heatmap is plotted, and a scatterplot is drawn if only 2 columns.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return A correlation matrix.
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
#' # Step 1: Normalize raw data
#' normalized_obj <- proteo.normalize(raw_data, metadata)
#'
#' # Step 2: Compute correlation matrix
#' corr_mat <- proteo.corr(normalized_obj, metadata)
#'
#' @references
#' Pearson correlation: Pearson 1895 <doi:10.1098/rspl.1895.0041>
#' Spearman correlation: Spearman 1904 <doi:10.1037/1082-989X.9.2.220>
#'
#' @export


proteo.corr <- function(normalized_data, metadata = NULL, group_level = TRUE,
                        method = "auto", plot = TRUE, help = FALSE) {

  if (help || missing(normalized_data)) {
    message("
Function proteo.corr()

Description:
  Computes pairwise correlations among samples or groups for normalized proteomics data.
  Automatically selects Pearson or Spearman correlation if 'auto', and provides heatmap
  and scatterplot visualizations.

Usage:
  proteo.corr(normalized_data, metadata = NULL, group_level = TRUE, method = 'auto', plot = TRUE)

Arguments:
  normalized_data  Proteomics object (proteoNorm) or data.frame with ProteinID column.
  metadata         Sample metadata (required if group_level = TRUE), with 'Sample' and 'Group'.
  group_level      If TRUE, correlations are computed at group level (averaging samples per group).
  method           Correlation method: 'pearson', 'spearman', 'kendall', or 'auto' (default).
  plot             Show heatmap and scatterplot? Default TRUE.
  help             Show this message and example usage.

Example (functional):
 raw_data <- data.frame(
  ProteinID = paste0('P', 1:5),
  Control1 = c(100, 200, 150, 300, 250),
  Control2 = c(110, 210, 160, 310, 260),
  Treatment1 = c(300, 100, 200, 150, 250),
  Treatment2 = c(310, 90, 210, 140, 260)
)

metadata <- data.frame(
  Sample = c('Control1'', 'Control2', 'Treatment1', 'Treatment2'),
  Group = c('Control', 'Control', 'Treatment', 'Treatment')
)

# Step 1: Normalize
normalized_obj <- proteo.normalize(raw_data, metadata)

# Step 2: Compute correlation matrix
corr_mat <- proteo.corr(normalized_obj, metadata)
")
    return(invisible(NULL))
  }

  if (group_level && is.null(metadata)) {
    stop("metadata must be provided if group_level = TRUE")
  }

  # --- Check required packages ---
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install 'pheatmap'")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'")

  # --- Extract matrix if proteoNorm object ---
  if ("proteoNorm" %in% class(normalized_data)) {
    df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
    normalized_data <- normalized_data[[df_names]]
  }

  if (!"ProteinID" %in% colnames(normalized_data)) stop("Missing 'ProteinID' column")

  expr_df <- dplyr::select(normalized_data, -ProteinID)
  rownames(expr_df) <- normalized_data$ProteinID

  # --- Aggregate by group if requested ---
  if (!is.null(metadata) && group_level) {
    if (!all(c("Sample","Group") %in% colnames(metadata))) stop("Metadata must contain 'Sample' and 'Group' columns")

    df_long <- tidyr::pivot_longer(normalized_data, -ProteinID, names_to = "Sample", values_to = "Abundance")
    df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

    df_grouped <- dplyr::summarise(
      dplyr::group_by(df_long, ProteinID, Group),
      Abundance = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

    df_grouped <- tidyr::pivot_wider(df_grouped, names_from = Group, values_from = Abundance)

    corr_df <- dplyr::select(df_grouped, -ProteinID)
  } else {
    corr_df <- expr_df
  }

  # --- Choose correlation method ---
  if (method == "auto") {
    normality <- apply(corr_df, 2, function(x) tryCatch(stats::shapiro.test(x)$p.value > 0.05, error = function(e) FALSE))
    has_ties <- any(apply(corr_df, 2, function(x) anyDuplicated(x) > 0))
    method_used <- if (all(normality) && !has_ties) "pearson" else "spearman"
    message("Automatically selected correlation method: ", method_used)
  } else {
    method_used <- method
  }

  # --- Correlation matrix ---
  corr_mat <- stats::cor(corr_df, method = method_used, use = "pairwise.complete.obs")

  # --- Create "r = ..." matrix for display ---
  corr_text <- matrix(paste0("r = ", round(corr_mat, 2)),
                      nrow = nrow(corr_mat),
                      ncol = ncol(corr_mat),
                      dimnames = dimnames(corr_mat))

  # --- Plot heatmap ---
  if (plot) {
    pheatmap::pheatmap(
      corr_mat,
      display_numbers = corr_text,
      number_color = "black",
      main = paste("Correlation (", method_used,")"),
      cluster_rows = FALSE,
      cluster_cols = FALSE
    )

    # Scatterplot if only 2 columns
    if (ncol(corr_df) == 2) {
      cols <- colnames(corr_df)
      test <- stats::cor.test(corr_df[[1]], corr_df[[2]], method = method_used)
      r_val <- round(test$estimate, 3)
      p_val <- if(test$p.value < 0.001) "<0.001" else signif(test$p.value, 3)

      strength <- cut(abs(r_val),
                      breaks = c(-Inf, 0.3, 0.5, 0.7, 0.9, Inf),
                      labels = c("very weak or none", "weak", "moderate", "strong", "very strong"),
                      right = FALSE)

      g <- ggplot2::ggplot(corr_df, ggplot2::aes(x = .data[[cols[1]]], y = .data[[cols[2]]])) +
        ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
        ggplot2::geom_smooth(method = ifelse(method_used=="pearson","lm","loess"),
                             se = FALSE, color = "red", linetype = "dashed") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Correlation scatterplot (", method_used, ")", sep=""),
                      subtitle = paste0("r = ", r_val, " | p = ", p_val, " | ", strength),
                      x = cols[1], y = cols[2])
      print(g)
    }
  }

  return(corr_mat)
}
