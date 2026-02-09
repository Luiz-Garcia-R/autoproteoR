#' Protein Ranking by Abundance or Variance
#'
#' Generates a per-group ranking of proteins based on their mean abundance or variance.
#' The function works with a `proteoNorm` object and can optionally highlight specific proteins in the plot.
#'
#' @param norm_data A `proteoNorm` object returned by `proteo.normalize()`.
#' @param group_col Character. Column name in `metadata` that defines sample groups (default: "Group").
#' @param metric Character. Metric to rank proteins by, either "mean" (default) or "var" (variance).
#' @param top_n Integer. If provided, only the top N proteins per group are retained.
#' @param highlight_protein Character vector. Names of proteins to highlight in the plot (default: NULL).
#' @param plot Logical. If TRUE (default), a ggplot2 ranking plot is displayed.
#'
#' @return Invisibly returns a list with two elements:
#'   \describe{
#'     \item{ranked_proteins}{Data frame with per-group protein statistics and ranks.}
#'     \item{plot}{ggplot2 object with protein rankings.}
#'   }
#'
#' @details
#' The function computes mean and variance for each protein in each group and assigns ranks.
#' It supports exactly two groups for comparison. Highlighted proteins are displayed with
#' larger points and labels.
#'
#' @examples
#' \dontrun{
#' @examples
#' # --- Small replicable dataset ---
#' normalized_matrix <- data.frame(
#'   ProteinID = paste0("P", 1:5),
#'   Control_1 = c(20, 22, 19, 25, 18),
#'   Control_2 = c(21, 23, 20, 24, 19),
#'   Treatment_1 = c(30, 28, 32, 27, 31),
#'   Treatment_2 = c(29, 27, 33, 26, 30)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control_1","Control_2","Treatment_1","Treatment_2"),
#'   Group  = c("Control","Control","Treatment","Treatment")
#' )
#'
#' proteo_obj <- list(normalized_matrix = normalized_matrix,
#'                    metadata = metadata)
#' class(proteo_obj) <- "proteoNorm"
#'
#' proteo.rank(proteo_obj, metric = "mean", top_n = 3, highlight_protein = "P1")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text facet_wrap labs theme_minimal
#' @export

proteo.rank <- function(normalized_data,
                        group_col = "Group",
                        metric = c("mean", "var"),
                        top_n = NULL,
                        highlight_protein = NULL,
                        plot = TRUE) {

  metric <- match.arg(metric)

  # --- Extract expression and metadata ---

  if ("proteonorm" %in% class(normalized_data)) {
    expr <- normalized_data$normalized_matrix
    meta <- normalized_data$metadata
  } else {
    stop("Input must be a 'proteonorm' object.")
  }

  # --- Check group column ---

  if (!group_col %in% colnames(meta)) stop("Group column not found in metadata.")
  groups <- unique(meta[[group_col]])
  if (length(groups) != 2) stop("Exactly 2 groups are required.")

  # --- Ensure ProteinID column exists ---

  if (!"ProteinID" %in% colnames(expr)) stop("Expression matrix must have a 'ProteinID' column")

  # --- Transform to long format for plotting ---

  library(dplyr)
  library(tidyr)
  expr_long <- expr %>%
    tidyr::pivot_longer(-ProteinID, names_to = "Sample", values_to = "Abundance") %>%
    dplyr::left_join(meta, by = "Sample")

  # --- Compute ranking metric per protein and group ---

  rank_df <- expr_long %>%
    dplyr::group_by(ProteinID, !!rlang::sym(group_col)) %>%
    dplyr::summarize(Metric = if(metric=="mean") mean(Abundance) else var(Abundance), .groups="drop") %>%
    dplyr::arrange(!!rlang::sym(group_col), dplyr::desc(Metric)) %>%
    dplyr::group_by(!!rlang::sym(group_col)) %>%
    dplyr::mutate(Rank = dplyr::row_number())

  # --- Filter top N if requested ---

  if(!is.null(top_n)){
    rank_df <- rank_df %>% dplyr::group_by(!!rlang::sym(group_col)) %>% dplyr::slice_head(n = top_n)
  }

  # --- Merge rank back to long data ---

  plot_df <- expr_long %>% dplyr::inner_join(rank_df %>% dplyr::select(ProteinID, !!rlang::sym(group_col), Rank),
                                             by = c("ProteinID", group_col))

  # --- Plot ---

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Rank, y = Abundance, color = !!rlang::sym(group_col))) +
    ggplot2::geom_point(alpha = 0.5, size = 2) +
    ggplot2::facet_wrap(as.formula(paste("~", group_col)), scales = "free_x", nrow = 1) +
    ggplot2::labs(x = "Protein rank", y = paste0("Abundance (", metric, ")"),
                  title = paste("Protein ranking by", metric, "per group")) +
    ggplot2::theme_minimal(base_size = 12)

  # --- Highlight specific proteins ---

  if(!is.null(highlight_protein)){
    highlight_points <- plot_df %>% dplyr::filter(ProteinID %in% highlight_protein)

    highlight_labels <- highlight_points %>% dplyr::group_by(ProteinID, !!rlang::sym(group_col)) %>%
      dplyr::slice_max(Abundance, n=1, with_ties = FALSE) %>% dplyr::ungroup()

    p <- p +
      ggplot2::geom_point(data=highlight_points, color="#FF6600", size=4) +
      ggplot2::geom_text(data=highlight_labels, ggplot2::aes(label=ProteinID),
                         vjust=-1, size=4, color="black")
  }

  if(plot) print(p)
  invisible(list(ranked_proteins = rank_df, plot = p))
}
