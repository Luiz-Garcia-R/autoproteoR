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
#'
#' @return Invisibly returns:
#'       \item{test}{The statistical test result (htest object).}
#'       \item{data}{Long-format data frame used for the test and plotting.}
#'       \item{plot}{A ggplot2 object with boxplot + jitter.}
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
#' # Perform limma moderated t-test for protein "P1"
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
                           protein_id = NULL,
                           group_col = "Group") {

  # -------------------------
  # Packages
  # -------------------------

  pkgs <- c("ggplot2", "dplyr", "limma", "tidyr")
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
  input_id <- protein_id

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
  rownames(mat) <- sub(".*\\|", "", rownames(mat))
  rownames(mat) <- sub("\\|.*", "", rownames(mat))

  mat <- mat[, !names(mat) %in% protein_col]
  mat <- as.matrix(mat)

  storage.mode(mat) <- "numeric"

  # -------------------------
  # Resolve protein identifier
  # -------------------------

  if (!(protein_id %in% rownames(mat))) {

    if (!requireNamespace("clusterProfiler", quietly = TRUE))
      stop("Package 'clusterProfiler' required for SYMBOL lookup.")

    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("Install 'org.Hs.eg.db' for SYMBOL mapping.")

    anno <- clusterProfiler::bitr(
      protein_id,
      fromType = "SYMBOL",
      toType = "UNIPROT",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )

    if (nrow(anno) == 0)
      stop("Protein not found as UNIPROT or SYMBOL.")

  # -------------------------
  # Map Uniprot ID
  # -------------------------
  mapped_ids <- unique(anno$UNIPROT)
  mapped_id <- intersect(mapped_ids, rownames(mat))

  if (length(mapped_id) == 0)
    stop("Mapped UNIPROT not present in dataset.")

  protein_id <- mapped_id[1]
  }

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

  groups <- levels(meta[[group_col]])

  g1 <- df_long$Abundance[df_long[[group_col]] == groups[1]]
  g2 <- df_long$Abundance[df_long[[group_col]] == groups[2]]

  plot_title <- protein_id
  if (exists("anno") && nrow(anno) > 0) {
    plot_title <- anno$SYMBOL[1]
  }

  # -------------------------
  # P label
  # -------------------------
  pval <- prot_res$P.Value

  signif_label <-
    if (pval < 0.001) "***"
  else if (pval < 0.01) "**"
  else if (pval < 0.05) "*"
  else "ns"

  p_label <- if (pval < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", formatC(pval, format = "f", digits = 3))
  }

  # -------------------------
  # Plot
  # -------------------------
  y_max <- max(df_long$Abundance, na.rm = TRUE)
  y_pos <- y_max * 1.08
  y_lim <- y_max * 1.1

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[group_col]], y = Abundance, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.5) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.65, color = "black") +
    ggplot2::annotate("segment",
                      x = 1, xend = 2,
                      y = y_pos * 0.98, yend = y_pos * 0.98) +
    ggplot2::annotate("text",
                      x = 1.5, y = y_pos,
                      label = signif_label,
                      size = 5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = paste0(
        "logFC = ", round(prot_res$logFC, 2),
        " | p = ", p_label),
      x = "",
      y = "Normalized abundance"
    ) +
    ggplot2::scale_y_continuous(limits = c(NA, y_lim), expand = c(0.02, 0))

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

  # --- Summary ---
  cat("=== Results ===\n")
  cat("Protein:", input_id, "\n")
  cat("Groups:", groups[1], "vs", groups[2], "\n")
  cat("Means:", signif(mean(g1, na.rm = TRUE), 3),
      "vs", signif(mean(g2, na.rm = TRUE), 3), "\n")
  cat("logFC:", signif(prot_res$logFC, 3), "\n")
  cat("p-value:", signif(pval, 3), signif_label, "\n")
  cat("================\n")

  return(invisible(out))

}
