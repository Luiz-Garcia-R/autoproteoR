#' Statistical test for protein abundance between two groups
#'
#' Performs a limma t moderated test for the normalized abundance of a specific protein across two groups.
#' Returns a significance summary and a boxplot with jittered points.
#'
#' @param normalized_data A normalized abundance data frame (containing `ProteinID`)
#'   or an object of class `proteonorm` (output of \code{proteo.normalize()}).
#' @param metadata A data frame with sample metadata. Must contain a `Sample` column
#'   and the grouping column specified in `group_col`.
#' @param protein_id Character. The ID of the protein to test (must be present in `normalized_data`).
#' @param group_col Character. Column name in `metadata` defining the groups
#'   (default = `"Group"`). Only two groups are currently supported.
#' @param return_type Either `"htest"` (returns only the test result) or `"all"`
#'   (returns test, data, and plot). Default = `"htest"`.
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
#' # Perform limma t-test for protein "P1"
#' result <- proteo.boxplot(normalized_data, protein_id = "P1")
#' }
#'
#' @references
#' Student, W. (1908). The probable error of a mean. \emph{Biometrika}, 6(1), 1–25. <doi:10.1093/biomet/6.1.1>
#' Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. \emph{Ann. Math. Stat.}, 18, 50–60. <doi:10.1214/aoms/1177730491>
#'
#' @export

proteo.boxplot <- function(normalized_data,
                           metadata = NULL,
                           protein_id,
                           group_col = "Group",
                           return_type = c("all", "htest")) {

  return_type <- match.arg(return_type)

  # -------------------------
  # Packages
  # -------------------------
  if (!requireNamespace("limma", quietly = TRUE))
    stop("Package 'limma' is required.")

  pkgs <- c("ggplot2", "dplyr", "tidyr")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Install: ", paste(miss, collapse = ", "))


  # -------------------------
  # proteonorm support
  # -------------------------
  if ("proteonorm" %in% class(normalized_data)) {

    if (is.null(normalized_data$normalized_matrix) ||
        is.null(normalized_data$metadata)) {

      stop("proteonorm must contain matrix + metadata")
    }

    metadata <- normalized_data$metadata
    normalized_data <- normalized_data$normalized_matrix
  }


  if (is.null(metadata))
    stop("metadata required")


  # -------------------------
  # Detect ProteinID
  # -------------------------
  id_names <- c("ProteinID", "Protein_ID", "Accession")

  protein_col <- intersect(id_names, names(normalized_data))[1]

  if (is.na(protein_col))
    stop("ProteinID column not found")


  # -------------------------
  # Sample column
  # -------------------------
  sample_names <- c("Sample", "sample")

  sample_col <- intersect(sample_names, names(metadata))[1]

  if (is.na(sample_col))
    stop("Sample column not found")


  # -------------------------
  # Group column
  # -------------------------
  if (!group_col %in% names(metadata))
    stop("group_col not found")


  # -------------------------
  # Build expression matrix
  # -------------------------
  mat <- normalized_data

  rownames(mat) <- mat[[protein_col]]

  mat <- mat[, !names(mat) %in% protein_col]

  mat <- as.matrix(mat)

  storage.mode(mat) <- "numeric"


  if (!protein_id %in% rownames(mat))
    stop("Protein not found")


  # -------------------------
  # Metadata order
  # -------------------------
  meta <- metadata

  meta <- meta[
    match(colnames(mat), meta[[sample_col]]),
  ]

  if (any(is.na(meta[[sample_col]])))
    stop("Sample mismatch")


  # -------------------------
  # Design (limma)
  # -------------------------
  meta[[group_col]] <- factor(meta[[group_col]])

  if (nlevels(meta[[group_col]]) != 2)
    stop("Exactly 2 groups required")


  design <- model.matrix(~ 0 + meta[[group_col]])
  colnames(design) <- levels(meta[[group_col]])

  contrast <- paste0(colnames(design)[2], "-", colnames(design)[1])

  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast,
    levels = design
  )


  # -------------------------
  # limma fit
  # -------------------------
  fit <- limma::lmFit(mat, design)

  fit <- limma::contrasts.fit(fit, contrast_matrix)

  fit <- limma::eBayes(fit)


  tt <- limma::topTable(
    fit,
    coef = 1,
    number = Inf,
    sort.by = "none"
  )


  prot_res <- tt[protein_id, ]


  # -------------------------
  # Long data for plot
  # -------------------------
  df_long <- data.frame(
    Sample = colnames(mat),
    Abundance = as.numeric(mat[protein_id, ])
  )

  df_long <- dplyr::left_join(
    df_long,
    meta,
    by = setNames(sample_col, "Sample")
  )


  # -------------------------
  # P label
  # -------------------------
  pval <- prot_res$P.Value

  signif_label <-
    if (pval < 0.001) "***"
  else if (pval < 0.01) "**"
  else if (pval < 0.05) "*"
  else "ns"


  y_pos <- max(df_long$Abundance, na.rm = TRUE) * 1.1


  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data[[group_col]],
      y = Abundance,
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.4, adjust = 0.6) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, color = "gray20", linewidth = 0.4) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1), alpha = 0.4, size = 1.8, color = "gray25") +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::annotate(
      "text", x = 1.5, y = y_pos,
      label = signif_label, size = 5
    ) +
    ggplot2::labs(
      title = protein_id,
      x = "",
      y = "Normalized abundance"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )


  print(p)

  # -------------------------
  # Output
  # -------------------------
  out <- list(
    logFC = prot_res$logFC,
    t     = prot_res$t,
    pval  = pval,
    adjp  = prot_res$adj.P.Val,
    plot  = p
  )


  if (return_type == "all") {
    return(out)
  } else {
    return(prot_res)
  }
}
