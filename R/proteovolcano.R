#' Volcano plot for proteomics analysis
#'
#' Generates a volcano plot from normalized proteomics data. Highlights
#' differentially expressed proteins (DEPs) based on adjusted p-value and
#' log2 fold-change thresholds. Protein labels can optionally be displayed
#' on the plot.
#'
#' @param normalized_data A data frame with normalized abundances
#'   (containing a `ProteinID` column) or an object of class `proteoNorm`.
#' @param group_col Column name in `metadata` that identifies the groups
#'   (default = `"Group"`). Only two groups are supported.
#' @param padj_threshold Adjusted p-value cutoff (default = 0.05).
#' @param log2fc_threshold Log2 fold-change cutoff (default = 1).
#' @param results Logical. If TRUE, returns a list with DEPs (default = FALSE).
#' @param identify Logical. If TRUE, adds protein labels directly on the plot (default = FALSE).
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

proteo.volcano <- function(normalized_data,
                           group_col = "Group",
                           padj_threshold = 0.05,
                           log2fc_threshold = 1,
                           results = FALSE,
                           identify = FALSE) {

  # --- Extract metadata and expression matrix if proteonorm object ---
  if ("proteonorm" %in% class(normalized_data)) {
    expr_df <- normalized_data$normalized_matrix
    metadata <- normalized_data$metadata
  } else {
    expr_df <- as.data.frame(normalized_data)
  }

  # --- Detect ProteinID column robustly ---
  possible_id_names <- c("ProteinID", "Protein_ID", "Protein", "Accession",
                         "Protein.Accession", "protein_id")
  protein_col <- intersect(names(expr_df), possible_id_names)[1]
  if (is.null(protein_col)) stop("Unable to detect ProteinID column in normalized_data.")
  protein_ids <- expr_df[[protein_col]]

  # --- Prepare expression matrix ---
  expr_cols <- setdiff(names(expr_df), protein_col)
  if (length(expr_cols) == 0) stop("No expression columns detected in normalized_data after removing ProteinID.")
  expr_mat <- as.matrix(expr_df[, expr_cols, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- protein_ids

  # --- Check metadata ---
  if (is.null(metadata)) stop("metadata must be provided or included in proteoNorm object.")
  possible_sample_names <- c("Sample", "sample", "SAMPLE")
  sample_col <- intersect(possible_sample_names, names(metadata))[1]
  if (is.null(sample_col)) stop("Metadata must contain a 'Sample' column.")

  # --- Group column detection ---
  if (is.null(group_col) || !group_col %in% names(metadata)) {
    possible_group_names <- c("Group", "group", "GROUP", "Condition", "condition")
    group_col <- intersect(possible_group_names, names(metadata))[1]
    if (is.null(group_col)) stop("Unable to detect Group column in metadata.")
  }

  # --- Extract samples per group ---
  groups <- unique(metadata[[group_col]])
  if (length(groups) != 2) stop("This function currently supports exactly 2 groups.")
  g1_samples <- metadata[[sample_col]][metadata[[group_col]] == groups[1]]
  g2_samples <- metadata[[sample_col]][metadata[[group_col]] == groups[2]]

  missing <- setdiff(c(g1_samples, g2_samples), colnames(expr_mat))
  if (length(missing) > 0) stop("Samples not found in normalized_data: ", paste(missing, collapse = ", "))

  expr_sub <- expr_mat[, c(g1_samples, g2_samples), drop = FALSE]

  # --- limma design ---
  group_factor <- factor(metadata[[group_col]][match(colnames(expr_sub), metadata[[sample_col]])])

  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)

  contrast <- paste0(levels(group_factor)[2], "-", levels(group_factor)[1])

  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast,
    levels = design
  )

  # --- limma fit ---
  fit <- limma::lmFit(expr_sub, design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit)

  tt <- limma::topTable(
    fit,
    coef = 1,
    number = Inf,
    sort.by = "none"
  )

  res_df <- data.frame(
    Protein = rownames(tt),
    log2FoldChange = tt$logFC,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val
  )

  # --- Annotate proteins via clusterProfiler ---
  # --- Annotate proteins (UNIPROT -> SYMBOL) ---

  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Package 'clusterProfiler' required for annotation.")

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Install 'org.Hs.eg.db' for human annotation.")

  clean_uniprot <- function(x) {
    x <- sub("^\\w+\\|", "", x, perl = TRUE)
    x <- sub("\\|.*$", "", x, perl = TRUE)
    x <- sub("-\\d+$", "", x, perl = TRUE)
    trimws(x)
  }

  res_df$UNIPROT <- clean_uniprot(res_df$Protein)

  anno <- clusterProfiler::bitr(
    unique(res_df$UNIPROT),
    fromType = "UNIPROT",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  res_df <- merge(
    res_df,
    anno[, c("UNIPROT", "SYMBOL")],
    by = "UNIPROT",
    all.x = TRUE
  )

  res_df$PlotLabel <- ifelse(
    !is.na(res_df$SYMBOL) & nzchar(res_df$SYMBOL),
    res_df$SYMBOL,
    res_df$UNIPROT
  )

  # --- Volcano summary ---
  up <- res_df$UNIPROT[
    res_df$padj < padj_threshold &
      res_df$log2FoldChange > log2fc_threshold
  ]

  down <- res_df$UNIPROT[
    res_df$padj < padj_threshold &
      res_df$log2FoldChange < -log2fc_threshold
  ]

  # --- Volcano plot ---
  res_df$group_color <- factor(ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange > log2fc_threshold, "Up",
                                      ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange < -log2fc_threshold, "Down", "NS")))

  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = group_color)) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("Up" = "#ff3333", "Down" = "#006699", "NS" = "darkgrey")) +
    ggplot2::labs(color = "Regulation",
                  x = "log2 Fold Change",
                  y = "-log10(p-adj)",
                  title = paste0(groups[2], " vs ", groups[1])) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))

  if (identify) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = PlotLabel),
      size = 3,
      max.overlaps = 50
    )
  }

  print(p)

  # --- Return ---
  message("=== Volcano summary ===")
  message("Groups: ", groups[1], " vs ", groups[2])
  message("Up (", groups[2], "): ", length(up))
  message("Down (", groups[1], "): ", length(down))
  message("padj threshold: ", padj_threshold, " | log2FC threshold: ", log2fc_threshold)
  message("=================================")

  volcano_deps <- list(up = up, down = down, full_results = res_df)

  if (results) return(volcano_deps)
  invisible(volcano_deps)
}
