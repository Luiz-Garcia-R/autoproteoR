#' Identify Differentially Expressed Proteins in Proteomics Data
#'
#' This function performs pairwise comparison between two groups of samples
#' in a normalized proteomics dataset. It calculates log2 fold-changes,
#' performs t-tests per protein, adjusts p-values for multiple testing,
#' and annotates UniProt IDs using \code{clusterProfiler::bitr}. Optionally,
#' it generates a bar plot of the top up- and down-regulated proteins.
#'
#' @param normalized_data A list returned by \code{proteo.normalize()} containing
#'   normalized protein expression in \code{x$normalized_matrix} and
#'   sample metadata in \code{x$metadata}.
#' @param group_col Character. Column in \code{x$metadata} specifying group labels. Default = "Group".
#' @param id_col Character. Column in \code{x$normalized_matrix} containing UniProt IDs. Default = "ProteinID".
#' @param p_adj_cutoff Numeric. Adjusted p-value threshold to define significance. Default = 0.05.
#' @param logfc_cutoff Numeric. Log2 fold-change threshold to define up- or down-regulation. Default = 0.5.
#' @param top_n Integer. Number of top up- and down-regulated proteins to display. Default = 10.
#' @param organism Character. Organism database for annotation: "hs" (human) or "mm" (mouse). Default = "hs".
#' @param plot Logical. If TRUE (default), generates a bar plot of top DEPs.
#' @param verbose Logical. If TRUE (default), prints messages and progress information.
#'
#' @return A list of class \code{proteoDeps} containing:
#' \describe{
#'   \item{summary}{Summary of comparison: number of proteins tested, number significant, counts of up/down-regulated.}
#'   \item{table}{Data frame with all proteins, log2 fold-change, p-value, adjusted p-value, and direction.}
#'   \item{significant}{Subset of \code{table} containing only significant DEPs.}
#'   \item{top}{Data frame with top up- and down-regulated proteins, cleaned UniProt IDs, mapped symbols, and plotting labels.}
#'   \item{plot}{ggplot2 object with top DEPs bar plot (or NULL if plot = FALSE).}
#'   \item{metadata}{List of function parameters and number of samples per group.}
#' }
#'
#' @examples
#' \dontrun{
#' raw_data <- data.frame(
#' ProteinID = c("P04637", "P31749", "Q16539", "P00533", "P38398"), # real human UniProt IDs
#' Control1 = c(100, 200, 150, 300, 250),
#' Control2 = c(110, 210, 160, 310, 260),
#' Treatment1 = c(300, 100, 200, 150, 250),
#' Treatment2 = c(310, 90, 210, 140, 260)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control1", "Control2", "Treatment1", "Treatment2"),
#'   Group = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' normalized_obj <- proteo.normalize(raw_data, metadata, method = "none")
#' proteo.deps(normalized_obj)
#' }
#'
#' @references
#' Benjamini & Hochberg, 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B 57:289-300.
#' UniProt: The Universal Protein Resource, 2021. Nucleic Acids Res 49:D394–D403.
#'
#' @export

proteo.deps <- function(normalized_data,
                        group_col = "Group",
                        id_col = "ProteinID",
                        p_adj_cutoff = 0.05,
                        logfc_cutoff = 0.5,
                        top_n = 10,
                        organism = c("mouse", "human"),
                        plot = TRUE,
                        verbose = TRUE) {

  organism <- match.arg(organism)

  # --- 1) Basic checks ---

  # Check required packages ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Please install 'clusterProfiler'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")

  if (is.null(normalized_data$normalized_matrix)) stop("Object does not contain `normalized_matrix`.")
  if (is.null(normalized_data$metadata)) stop("Object does not contain `metadata`.")
  if (!group_col %in% colnames(normalized_data$metadata)) stop(paste("Group column not found:", group_col))
  if (!id_col %in% colnames(normalized_data$normalized_matrix)) stop(paste("ID column not found:", id_col))

  mat <- as.data.frame(normalized_data$normalized_matrix, stringsAsFactors = FALSE)
  meta <- as.data.frame(normalized_data$metadata, stringsAsFactors = FALSE)

  # Ensure sample names match

  sample_cols <- intersect(colnames(mat), meta$Sample)
  if (length(sample_cols) == 0) stop("No matching sample names found between matrix and metadata$Sample.")

  # --- 2) Define groups automatically ---

  groups <- unique(meta[[group_col]])
  if (length(groups) != 2) stop("Exactly two groups are required.")

  g1 <- groups[1]; g2 <- groups[2]
  if (verbose) message("Comparing groups: ", g1, " vs ", g2)

  samples_g1 <- meta$Sample[meta[[group_col]] == g1]
  samples_g2 <- meta$Sample[meta[[group_col]] == g2]

  # --- 3) T-test per protein ---

  do_ttest <- function(i) {
    x1 <- as.numeric(mat[i, samples_g1])
    x2 <- as.numeric(mat[i, samples_g2])
    if (all(is.na(x1)) || all(is.na(x2))) return(c(logFC = NA, pvalue = NA))
    logFC <- mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE)
    pval <- tryCatch(stats::t.test(x1, x2)$p.value, error = function(e) NA_real_)
    c(logFC = logFC, pvalue = pval)
  }

  stats <- t(vapply(seq_len(nrow(mat)), do_ttest, numeric(2)))
  stats <- as.data.frame(stats, stringsAsFactors = FALSE)
  stats[[id_col]] <- mat[[id_col]]  # add UniProt IDs
  stats$pvalue <- as.numeric(stats$pvalue)
  stats$padj <- p.adjust(stats$pvalue, method = "BH")

  # --- 4) Filter DEPs ---

  deps <- stats[!is.na(stats$logFC) & !is.na(stats$padj), ]
  deps$direction <- ifelse(deps$padj < p_adj_cutoff & deps$logFC >= logfc_cutoff, "Up",
                           ifelse(deps$padj < p_adj_cutoff & deps$logFC <= -logfc_cutoff, "Down", "NS"))
  deps_sig <- deps[deps$direction != "NS", ]

  if (verbose) message("Significant DEPs: ", nrow(deps_sig))

  if (nrow(deps_sig) == 0) {
    warning("No significant DEPs found.")
    return(list(summary = list(n_tested = nrow(deps), n_significant = 0),
                table = deps,
                significant = deps_sig,
                plot = NULL,
                metadata = list(p_adj_cutoff = p_adj_cutoff, logfc_cutoff = logfc_cutoff,
                                groups = c(g1, g2), n_samples = table(meta[[group_col]]))))
  }

  # --- 5) Select top N up/down ---

  top_up <- head(deps_sig[deps_sig$direction == "Up", ][order(deps_sig$padj[deps_sig$direction == "Up"]), ], top_n)
  top_down <- head(deps_sig[deps_sig$direction == "Down", ][order(deps_sig$padj[deps_sig$direction == "Down"]), ], top_n)
  top <- rbind(top_up, top_down)

  # --- 6) Annotate proteins via clusterProfiler ---

  OrgDb <- switch(organism, mouse = org.Mm.eg.db::org.Mm.eg.db, human = org.Hs.eg.db::org.Hs.eg.db)

  clean_uniprot <- function(x) {
    x <- sub("^\\w+\\|", "", x, perl = TRUE)
    x <- sub("\\|.*$", "", x, perl = TRUE)
    x <- sub("-\\d+$", "", x, perl = TRUE)
    trimws(x)
  }

  top$Protein_clean <- clean_uniprot(top[[id_col]])

  if (verbose) {
    message("IDs for mapping via bitr():")
    print(head(top$Protein_clean, 20))
  }

  anno <- clusterProfiler::bitr(
    top$Protein_clean,
    fromType = "UNIPROT",
    toType = "SYMBOL",
    OrgDb = OrgDb
  )

  top <- merge(top, anno[, c("UNIPROT", "SYMBOL")],
               by.x = "Protein_clean", by.y = "UNIPROT", all.x = TRUE)

  top$PlotLabel <- ifelse(!is.na(top$SYMBOL) & nzchar(top$SYMBOL), top$SYMBOL, top[[id_col]])

  # --- 7) Prepare dataframe for plotting ---

  top$PlotLabel <- factor(top$PlotLabel, levels = top$PlotLabel[order(top$logFC)])

  # --- 8) Plot ---

  p <- NULL
  if (plot) {
    p <- ggplot2::ggplot(top, ggplot2::aes(y = PlotLabel, x = logFC, fill = direction)) +
      ggplot2::geom_col() +
      ggplot2::labs(y = "", x = "log2 fold-change",
                    title = sprintf("Top %d Up- and Down-regulated Proteins (%s vs %s)", top_n, g2, g1)) +
      ggplot2::scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    print(p)
  }

  # --- 9) Summary and return ---

  summary_list <- list(
    comparison = paste0(g1, " vs ", g2),
    n_tested = nrow(deps),
    n_significant = nrow(deps_sig),
    n_up = sum(deps_sig$direction == "Up"),
    n_down = sum(deps_sig$direction == "Down"),
    top = top
  )

  output <- list(
    summary = summary_list,
    table = deps,
    significant = deps_sig,
    top = top,
    plot = p,
    metadata = list(p_adj_cutoff = p_adj_cutoff,
                    logfc_cutoff = logfc_cutoff,
                    groups = c(g1, g2),
                    n_samples = table(meta[[group_col]]))
  )

  class(output) <- "proteoDeps"
  invisible(output)
}
