#' Volcano plot for proteomics analysis
#'
#' Generates a volcano plot from normalized proteomics data. Highlights
#' differentially expressed proteins (DEPs) based on adjusted p-value and
#' log2 fold-change thresholds. Protein labels can optionally be displayed
#' on the plot.
#'
#' @param normalized_data A data frame with normalized abundances
#'   (containing a `ProteinID` column) or an object of class `proteoNorm`.
#' @param metadata A data frame with sample metadata. Must include columns
#'   `Sample` and the one defined in `group_col`.
#' @param group_col Column name in `metadata` that identifies the groups
#'   (default = `"Group"`). Only two groups are supported.
#' @param padj_threshold Adjusted p-value cutoff (default = 0.05).
#' @param log2fc_threshold Log2 fold-change cutoff (default = 1).
#' @param results Logical. If TRUE, returns a list with DEPs (default = FALSE).
#' @param identify Logical. If TRUE, adds protein labels directly on the plot (default = FALSE).
#' @param help Logical. If TRUE, prints this help message instead of running the function.
#'
#' @return Invisibly returns a list with three elements:
#'   \describe{
#'     \item{up}{Proteins classified as up-regulated.}
#'     \item{down}{Proteins classified as down-regulated.}
#'     \item{full_results}{Data frame with complete test results.}
#'   }
#'
#' @examples
#' \dontrun{
#' # --- Small example dataset ---
#' normalized_data <- data.frame(
#'   ProteinID = paste0("P", 1:6),
#'   Control1 = c(100, 200, 150, 80, 120, 160),
#'   Control2 = c(110, 210, 140, 90, 130, 170),
#'   Treatment1 = c(300, 100, 200, 60, 220, 180),
#'   Treatment2 = c(310, 90, 210, 70, 230, 190)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control1", "Control2", "Treatment1", "Treatment2"),
#'   Group  = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' # Step 1: Normalize (optional if already normalized)
#' normalized_obj <- proteo.normalize(normalized_data, metadata)
#'
#' # Step 2: Generate volcano plot
#' proteo.volcano(normalized_obj, metadata,
#'                group_col = "Group",
#'                padj_threshold = 0.05,
#'                log2fc_threshold = 1,
#'                results = TRUE,
#'                identify = TRUE)
#' }
#'
#' @references
#' Cui, X. & Churchill, G.A. (2003). Statistical tests for differential expression in
#'   cDNA microarray experiments. \emph{Genome Biology}, 4(4), 210. <doi:10.1186/gb-2003-4-4-210>
#' Khatri, P., Sirota, M., & Butte, A.J. (2012). Ten years of pathway analysis: current approaches
#'   and outstanding challenges. \emph{PLoS Comput Biol}, 8(2), e1002375. <doi:10.1371/journal.pcbi.1002375>
#'
#' @importFrom stats t.test p.adjust
#' @importFrom ggrepel geom_text_repel
#' @export

proteo.volcano <- function(normalized_data, metadata, group_col = "Group",
                           padj_threshold = 0.05, log2fc_threshold = 1,
                           results = FALSE, identify = FALSE, help = FALSE) {

  # --- Help message ---
  if (help || missing(normalized_data) || missing(metadata)) {
    message("
Function proteo.volcano()

Description:
  Generates a volcano plot from normalized proteomics data.
  Highlights differentially expressed proteins (DEPs) based on adjusted p-value
  and log2 fold-change thresholds. Protein labels can optionally be displayed.

Usage:
  proteo.volcano(normalized_data, metadata, group_col = 'Group',
                 padj_threshold = 0.05, log2fc_threshold = 1,
                 results = FALSE, identify = FALSE)

Arguments:
  normalized_data  Data frame with normalized abundances or proteoNorm object.
  metadata         Data frame with sample metadata. Must contain 'Sample' and group column.
  group_col        Column name in metadata identifying the groups. Default: 'Group'.
  padj_threshold   Adjusted p-value cutoff. Default: 0.05.
  log2fc_threshold Log2 fold-change cutoff. Default: 1.
  results          Logical. Return list with DEPs? Default: FALSE.
  identify         Logical. Add protein labels to plot? Default: FALSE.
  help             Logical. If TRUE, prints this help message.

Return:
  Invisibly returns a list with 'up', 'down', and 'full_results'.

Example:
  normalized_data <- data.frame(
    ProteinID = paste0('P',1:6),
    Control1 = c(100,200,150,80,120,160),
    Control2 = c(110,210,140,90,130,170),
    Treatment1 = c(300,100,200,60,220,180),
    Treatment2 = c(310,90,210,70,230,190)
  )
  metadata <- data.frame(
    Sample = c('Control1','Control2','Treatment1','Treatment2'),
    Group  = c('Control','Control','Treatment','Treatment')
  )
  proteo.volcano(normalized_data, metadata, group_col='Group', results=TRUE, identify=TRUE)
")
    return(invisible(NULL))
  }

  # --- Detect proteoNorm object ---
  if ("proteoNorm" %in% class(normalized_data)) {
    df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
    if (length(df_names) != 1) stop("Could not detect the normalized dataframe in the proteoNorm object.")
    normalized_data <- normalized_data[[df_names]]
  }

  # --- Prepare expression matrix ---
  if (!"ProteinID" %in% colnames(normalized_data)) {
    stop("normalized_data must contain a 'ProteinID' column.")
  }
  protein_ids <- normalized_data$ProteinID
  expr_mat <- normalized_data[, -which(colnames(normalized_data) == "ProteinID")]
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- protein_ids

  # --- Prepare metadata ---
  if (!"Sample" %in% colnames(metadata)) stop("Metadata must contain the column 'Sample'.")
  if (!group_col %in% colnames(metadata)) stop("Group column not found in metadata.")

  groups <- unique(metadata[[group_col]])
  if (length(groups) != 2) stop("This function currently supports exactly 2 groups.")

  g1 <- groups[1]
  g2 <- groups[2]

  g1_samples <- metadata$Sample[metadata[[group_col]] == g1]
  g2_samples <- metadata$Sample[metadata[[group_col]] == g2]

  missing <- setdiff(c(g1_samples, g2_samples), colnames(expr_mat))
  if (length(missing) > 0) {
    stop("Samples not found in normalized_data: ", paste(missing, collapse = ", "))
  }

  expr_sub <- expr_mat[, c(g1_samples, g2_samples), drop = FALSE]

  # --- Compute t-tests and log2 fold changes ---
  res_df <- data.frame(Protein = rownames(expr_sub),
                       log2FoldChange = NA_real_,
                       pvalue = NA_real_)

  for (i in seq_len(nrow(expr_sub))) {
    vals1 <- as.numeric(expr_sub[i, g1_samples])
    vals2 <- as.numeric(expr_sub[i, g2_samples])
    ttest <- try(stats::t.test(vals2, vals1), silent = TRUE)
    if (inherits(ttest, "try-error")) {
      res_df$pvalue[i] <- NA
      res_df$log2FoldChange[i] <- NA
    } else {
      res_df$pvalue[i] <- ttest$p.value
      res_df$log2FoldChange[i] <- mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE)
    }
  }

  res_df$padj <- stats::p.adjust(res_df$pvalue, method = "BH")
  res_df <- dplyr::filter(res_df, !is.na(pvalue))

  # --- Volcano summary ---
  up <- res_df$Protein[res_df$padj < padj_threshold & res_df$log2FoldChange > log2fc_threshold]
  down <- res_df$Protein[res_df$padj < padj_threshold & res_df$log2FoldChange < -log2fc_threshold]

  message("=== Volcano summary ===")
  message("Groups: ", g1, " vs ", g2)
  message("Up (", g2, "): ", length(up))
  message("Down (", g1, "): ", length(down))
  message("padj threshold: ", padj_threshold, " | log2FC threshold: ", log2fc_threshold)
  message("=================================")

  if (!results) {
    message("Tip: To see differentially expressed proteins, call with 'results = TRUE'.")
    message("You can also label proteins on the volcano plot using 'identify = TRUE'.")
  }

  # --- Volcano plot ---
  res_df$group_color <- factor(ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange > log2fc_threshold, "Up",
                                      ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange < -log2fc_threshold, "Down", "NS")))

  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue), color = group_color)) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("Up" = "#ff3333", "Down" = "#006699", "NS" = "darkgrey")) +
    ggplot2::labs(color = "Regulation",
                  x = "log2 Fold Change",
                  y = "-log10(p-value)",
                  title = paste0(g2, " vs ", g1)) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))

  # --- Add labels if identify = TRUE ---
  if (identify) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Protein), size = 3, max.overlaps = 50, na.rm = TRUE)
  }

  print(p)

  # --- Return DEPs if requested ---
  volcano_deps <- list(up = up, down = down, full_results = res_df)

  if (results) return(volcano_deps)
  invisible(volcano_deps)
}














