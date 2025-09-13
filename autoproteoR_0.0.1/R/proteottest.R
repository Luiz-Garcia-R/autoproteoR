#' Statistical test for protein abundance between two groups
#'
#' Performs a Student's t-test (default) or a Mann-Whitney U test (optional)
#' for the normalized abundance of a specific protein across two groups.
#' Returns a significance summary and a boxplot with jittered points.
#'
#' @param normalized_data A normalized abundance data frame (containing `ProteinID`)
#'   or an object of class `proteoNorm` (output of \code{proteo.normalize()}).
#' @param metadata A data frame with sample metadata. Must contain a `Sample` column
#'   and the grouping column specified in `group_col`.
#' @param protein_id Character. The ID of the protein to test (must be present in `normalized_data`).
#' @param group_col Character. Column name in `metadata` defining the groups
#'   (default = `"Group"`). Only two groups are currently supported.
#' @param utest Logical. If `TRUE`, performs Mann-Whitney test instead of t-test
#'   (default = `FALSE`). Use with caution on normalized data.
#' @param return_type Either `"htest"` (returns only the test result) or `"all"`
#'   (returns test, data, and plot). Default = `"htest"`.
#' @param help Logical. If TRUE, prints this help message instead of running the function.
#'
#' @return Invisibly returns:
#'   - if `return_type = "htest"`: an object of class \code{htest} (t-test or Wilcoxon test).
#'   - if `return_type = "all"`: a list with:
#'     \describe{
#'       \item{test}{The statistical test result (htest object).}
#'       \item{data}{Long-format data frame used for the test and plotting.}
#'       \item{plot}{A ggplot2 object with boxplot + jitter.}
#'     }
#'
#' @examples
#' \dontrun{
#' # Small example: normalized data
#' normalized_data <- data.frame(
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
#' # Perform t-test for protein "P1"
#' result <- proteo.ttest(normalized_data, metadata, protein_id = "P1", return_type = "all")
#' }
#'
#' @references
#' Student, W. (1908). The probable error of a mean. \emph{Biometrika}, 6(1), 1–25. <doi:10.1093/biomet/6.1.1>
#' Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. \emph{Ann. Math. Stat.}, 18, 50–60. <doi:10.1214/aoms/1177730491>
#'
#' @export

proteo.ttest <- function(normalized_data, metadata, protein_id, group_col = "Group",
                         utest = FALSE, return_type = c("htest", "all"), help = FALSE) {

  # --- Help message ---
  if (help || missing(normalized_data) || missing(metadata) || missing(protein_id)) {
    message("
Function proteo.ttest()

Description:
  Performs a t-test (default) or Mann-Whitney U test (optional) for a single protein
  across two groups in normalized proteomics data. Returns test results and optional boxplot.

Usage:
  proteo.ttest(normalized_data, metadata, protein_id, group_col='Group', utest=FALSE, return_type='htest')

Arguments:
  normalized_data  Normalized data frame (ProteinID column) or proteoNorm object.
  metadata         Data frame with sample metadata. Must contain 'Sample' and group_col.
  protein_id       Protein ID to test (must exist in normalized_data).
  group_col        Column in metadata defining groups (default 'Group'). Only 2 groups supported.
  utest            Logical. If TRUE, performs Mann-Whitney test instead of t-test. Default FALSE.
  return_type      'htest' returns only test result; 'all' returns list with test, data, plot.
  help             Logical. If TRUE, prints this help message.

Return:
  Invisibly returns:
    - if return_type='htest': htest object
    - if return_type='all': list with test, data, and ggplot2 boxplot

Example:
  normalized_data <- data.frame(
    ProteinID = paste0('P', 1:5),
    Control1 = c(100,200,150,300,250),
    Control2 = c(110,210,160,310,260),
    Treatment1 = c(300,100,200,150,250),
    Treatment2 = c(310,90,210,140,260)
  )
  metadata <- data.frame(
    Sample = c('Control1','Control2','Treatment1','Treatment2'),
    Group = c('Control','Control','Treatment','Treatment')
  )
  result <- proteo.ttest(normalized_data, metadata, protein_id='P1', return_type='all')
")
    return(invisible(NULL))
  }

  # --- Required packages ---
  pkgs <- c("ggplot2", "dplyr", "tidyr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Please install the packages: ", paste(missing_pkgs, collapse = ", "),
         ". Ex.: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))")
  }

  return_type <- match.arg(return_type)

  # --- Extract normalized matrix if proteoNorm object ---
  if ("proteoNorm" %in% class(normalized_data)) {
    df_names <- grep("_normalized$", names(normalized_data), value = TRUE)
    if (length(df_names) != 1) stop("Could not detect the normalized data.frame in 'proteoNorm'.")
    normalized_data <- normalized_data[[df_names]]
  }

  # --- Basic checks ---
  if (!"ProteinID" %in% colnames(normalized_data)) stop("Data.frame must contain 'ProteinID'.")
  if (!"Sample" %in% colnames(metadata)) stop("'metadata' must contain 'Sample'.")
  if (!group_col %in% colnames(metadata)) stop(paste0("Group column '", group_col, "' not found in metadata."))
  if (!protein_id %in% normalized_data$ProteinID) stop(paste0("Protein '", protein_id, "' not found in normalized_data."))

  # --- Extract protein values and join with metadata ---
  df <- normalized_data[normalized_data$ProteinID == protein_id, , drop = FALSE]
  df <- df[, setdiff(colnames(df), "ProteinID"), drop = FALSE]
  df_long <- tidyr::pivot_longer(df, cols = everything(), names_to = "Sample", values_to = "Abundance")
  df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

  # --- Ensure exactly 2 groups ---
  groups <- unique(df_long[[group_col]])
  if (length(groups) != 2) stop(paste("Test requires exactly 2 groups (found:", length(groups), ")."))

  # --- Numeric vectors per group ---
  g1 <- df_long$Abundance[df_long[[group_col]] == groups[1]]
  g1 <- g1[is.finite(g1)]

  g2 <- df_long$Abundance[df_long[[group_col]] == groups[2]]
  g2 <- g2[is.finite(g2)]

  if (length(g1) < 2 || length(g2) < 2) stop("Each group must have at least 2 non-NA values.")

  # --- Check for many zeros ---
  pct_zero_g1 <- mean(g1 == 0, na.rm = TRUE)
  pct_zero_g2 <- mean(g2 == 0, na.rm = TRUE)
  if (pct_zero_g1 > 0.5 || pct_zero_g2 > 0.5) {
    stop("More than 50% of values in one of the groups are zero. t-test or Mann-Whitney is not recommended.")
  }

  # --- Statistical test ---
  if (utest) {
    warning("You are running Mann-Whitney on normalized data. This is not recommended.")
    result <- stats::wilcox.test(g1, g2)
    test_name <- "Mann-Whitney"
  } else {
    result <- stats::t.test(g1, g2)
    test_name <- "t-test"
  }

  # --- Significance labels ---
  pval <- result$p.value
  signif_label <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"
  p_label <- if (pval < 0.001) "p < 0.001" else paste0("p = ", signif(pval, 3))

  # --- Plot ---
  y_max <- max(df_long$Abundance, na.rm = TRUE)
  y_pos <- y_max * 1.2
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[group_col]], y = Abundance, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.65, color = "black") +
    ggplot2::annotate("text", x = 1.5, y = y_pos, label = signif_label, size = 6) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(title = paste(test_name, "-", protein_id),
                  subtitle = p_label, x = "", y = "Normalized abundance") +
    ggplot2::scale_y_continuous(expand = c(0.02, 0.1))

  print(p)

  # --- Package-style message ---
  cat("=== Results (", test_name, ") ===\n", sep = "")
  cat("Protein:", protein_id, "\n")
  cat("Groups:", groups[1], "vs", groups[2], "\n")
  cat("Means:", signif(mean(g1), 3), "vs", signif(mean(g2), 3), "\n")
  cat("p-value:", signif(pval, 3), signif_label, "\n")
  cat("================\n")

  if (return_type == "all") invisible(list(test = result, data = df_long, plot = p))
  else invisible(result)
}






