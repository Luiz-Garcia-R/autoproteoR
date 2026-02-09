#' Statistical test for protein abundance between two groups
#'
#' Performs a Student's t-test (default) or a Mann-Whitney U test (optional)
#' for the normalized abundance of a specific protein across two groups.
#' Returns a significance summary and a boxplot with jittered points.
#'
#' @param normalized_data A normalized abundance data frame (containing `ProteinID`)
#'   or an object of class `proteonorm` (output of \code{proteo.normalize()}).
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

proteo.ttest <- function(normalized_data, metadata = NULL, protein_id, group_col = "Group",
                         utest = FALSE, return_type = c("htest", "all"), help = FALSE) {

  return_type <- match.arg(return_type)

  # --- Required packages ---
  pkgs <- c("ggplot2", "dplyr", "tidyr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse = ", "))

  return_type <- match.arg(return_type)

  # --- Extract metadata automatically if proteonorm object ---
  if ("proteonorm" %in% class(normalized_data)) {
    if (is.null(normalized_data$normalized_matrix) || is.null(normalized_data$metadata))
      stop("proteonorm object must contain $normalized_matrix and $metadata.")
    metadata <- normalized_data$metadata
    normalized_data <- normalized_data$normalized_matrix
  }

  # --- check metadata ---
  if (is.null(metadata)) stop("metadata must be provided")

  # --- Detect ProteinID column ---
  possible_id_names <- c("ProteinID", "Protein_ID", "Protein", "Accession", "Protein.Accession", "protein_id")
  protein_col <- names(normalized_data)[names(normalized_data) %in% possible_id_names][1]
  if (is.null(protein_col)) stop("Unable to detect ProteinID column.")

  # --- Detect Sample column ---
  if (is.null(metadata)) stop("metadata must be provided")
  possible_sample_names <- c("Sample", "sample", "SAMPLE")
  sample_col <- intersect(possible_sample_names, names(metadata))
  if (length(sample_col) == 0) stop("Unable to detect Sample column in metadata.")
  sample_col <- sample_col[1]

  # --- Use provided group_col if exists, else auto-detect ---
  if (is.null(group_col) || !group_col %in% names(metadata)) {
    possible_group_names <- c("Group", "group", "GROUP", "Condition", "condition")
    group_col <- intersect(possible_group_names, names(metadata))
    if (length(group_col) == 0) stop("Unable to detect Group column in metadata.")
    group_col <- group_col[1]
  }

  # --- Basic checks (mover PARA DEPOIS da extração do data.frame) ---
  if (!"ProteinID" %in% colnames(normalized_data)) stop("Data.frame must contain 'ProteinID'.")
  if (!"Sample" %in% colnames(metadata)) stop("'metadata' must contain 'Sample'.")
  if (!group_col %in% colnames(metadata)) stop(paste0("Group column '", group_col, "' not found in metadata."))
  if (!protein_id %in% normalized_data$ProteinID) stop(paste0("Protein '", protein_id, "' not found in normalized_data."))

  # --- Extract protein values ---
  df <- normalized_data[normalized_data[[protein_col]] == protein_id, , drop = FALSE]
  df <- df[, setdiff(names(df), protein_col), drop = FALSE]
  df_long <- tidyr::pivot_longer(df, cols = everything(), names_to = "Sample", values_to = "Abundance")
  df_long <- dplyr::left_join(df_long, metadata, by = setNames(sample_col, "Sample"))

  # --- Ensure exactly 2 groups ---
  groups <- unique(df_long[[group_col]])
  if (length(groups) != 2) stop(paste("Test requires exactly 2 groups (found:", length(groups), ")."))

  # --- Numeric vectors per group ---
  g1 <- df_long$Abundance[df_long[[group_col]] == groups[1]] |> na.omit()
  g2 <- df_long$Abundance[df_long[[group_col]] == groups[2]] |> na.omit()
  if (length(g1) < 2 || length(g2) < 2) stop("Each group must have at least 2 non-NA values.")

  # --- Check for many zeros ---
  if (mean(g1 == 0, na.rm = TRUE) > 0.5 || mean(g2 == 0, na.rm = TRUE) > 0.5) {
    stop("More than 50% of values in one of the groups are zero. Test not recommended.")
  }

  # --- Statistical test ---
  if (utest) {
    warning("You are running Mann-Whitney on normalized data. Not recommended.")
    result <- stats::wilcox.test(g1, g2)
    test_name <- "Mann-Whitney"
  } else {
    result <- stats::t.test(g1, g2)
    test_name <- "t-test"
  }

  # --- Plot ---
  y_max <- max(df_long$Abundance, na.rm = TRUE)
  y_pos <- y_max * 1.2
  pval <- result$p.value
  signif_label <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"
  p_label <- if (pval < 0.001) "p < 0.001" else paste0("p-value = ", signif(pval, 3))

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[group_col]], y = Abundance, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.65, color = "black") +
    ggplot2::annotate("text", x = 1.5, y = y_pos, label = signif_label, size = 5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(title = paste(protein_id),
                  subtitle = paste0(test_name, ": ", p_label), x = "", y = "Normalized abundance") +
    ggplot2::scale_y_continuous(expand = c(0.02, 0.1))

  print(p)

  # --- Summary ---
  cat("=== Results (", test_name, ") ===\n", sep = "")
  cat("Protein:", protein_id, "\n")
  cat("Groups:", groups[1], "vs", groups[2], "\n")
  cat("Means:", signif(mean(g1), 3), "vs", signif(mean(g2), 3), "\n")
  cat("p-value:", signif(pval, 3), signif_label, "\n")
  cat("================\n")

  if (return_type == "all") invisible(list(test = result, data = df_long, plot = p))
  else invisible(result)
}
