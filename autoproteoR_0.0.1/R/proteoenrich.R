#' Perform GO enrichment analysis for proteomics data
#'
#' Performs Gene Ontology (GO) enrichment analysis on proteomics data,
#' using either all detected proteins per group or only the unique proteins.
#' Supports human and mouse organisms and returns enriched GO terms with optional plots.
#'
#' @param normalized_data A proteomics object containing a `presence_per_group` element
#'   or a data.frame with a `ProteinID` column.
#' @param ont Character. Ontology category for GO analysis: `"BP"` (Biological Process),
#'   `"MF"` (Molecular Function), or `"CC"` (Cellular Component). Default = `"BP"`.
#' @param results Logical. If TRUE, returns a list with enrichment results and proteins used. Default = FALSE.
#' @param min_genes Minimum number of proteins required in a group to perform enrichment. Default = 3.
#' @param db Character. Organism database: `"hs"` (human, default) or `"mm"` (mouse).
#' @param show_plot Logical. If TRUE (default), shows dotplots for each group with significant GO terms.
#' @param p_adjust_cutoff Numeric. FDR threshold for filtering GO terms. Default = 0.05.
#' @param proteins Character. Which proteins to use per group: `"all"` (default) or `"unique"`.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#' @param verbose Logical. If TRUE (default), prints messages about progress.
#'
#' @return
#' If `results = TRUE`, returns a list with:
#' \describe{
#'   \item{ego_list}{List of `enrichResult` objects per group.}
#'   \item{proteins_used}{List of proteins used for enrichment per group.}
#' }
#' Otherwise, prints messages and optionally plots dotplots. Invisibly returns `NULL`.
#'
#' @examples
#' \dontrun{
#' # Small example: real human protein IDs
#' raw_data <- data.frame(
#'   ProteinID = c("P31946", "P62258", "P63104", "Q9Y262", "P62993"),
#'   Control_1  = c(1,1,0,1,0),
#'   Control_2  = c(1,1,0,1,0),
#'   Treatment_1 = c(0,1,1,1,1),
#'   Treatment_2 = c(0,1,1,1,1)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c('Control_1','Control_2','Treatment_1','Treatment_2'),
#'   Group  = c('Control','Control','Treatment','Treatment')
#' )
#'
#' # Step 1: Normalize data
#' normalized_obj <- proteo.normalize(raw_data, metadata)
#'
#' # Step 2: Run enrichment (min_genes = 1 for small example)
#' proteo.enrich(normalized_obj, min_genes = 1)
#' res <- proteo.enrich(normalized_obj, results = TRUE, min_genes = 1)
#' }
#'
#' @references
#' Yu et al. 2012. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16:284-287. <doi:10.1089/omi.2011.0118>
#'
#' @export

proteo.enrich <- function(normalized_data,
                          ont = "BP",
                          results = FALSE,
                          min_genes = 3,
                          db = c("hs", "mm"),
                          show_plot = TRUE,
                          p_adjust_cutoff = 0.05,
                          proteins = "all",
                          help = FALSE,
                          verbose = TRUE) {

  if (help || missing(normalized_data)) {
    message("
Function proteo.enrich()

Description:
  Performs Gene Ontology (GO) enrichment on proteomics data using all detected
  proteins or only unique proteins per group. Supports human and mouse organisms
  and returns enrichment results with optional dotplots.

Usage:
  proteo.enrich(normalized_data, ont = 'BP', ...)

Arguments:
  normalized_data  Proteomics object with 'presence_per_group' or data.frame with ProteinID.
  ont              GO category: 'BP', 'MF', or 'CC'. Default 'BP'.
  results          Return results as a list? Default FALSE.
  min_genes        Minimum proteins per group to run enrichment. Default 3.
  db               Organism: 'hs' (human) or 'mm' (mouse). Default 'hs'.
  show_plot        Show dotplots? Default TRUE.
  p_adjust_cutoff  FDR threshold for significant terms. Default 0.05.
  proteins         'all' or 'unique' proteins per group. Default 'all'.
  help             Show this message.

# Example dataset with real human protein IDs
raw_data <- data.frame(
  ProteinID = c('P31946', 'P62258', 'P63104', 'Q9Y262', 'P62993'), # real UniProt IDs
  Control_1  = c(1, 1, 0, 1, 0),
  Control_2  = c(1, 1, 0, 1, 0),
  Treatment_1 = c(0, 1, 1, 1, 1),
  Treatment_2 = c(0, 1, 1, 1, 1)
)

metadata <- data.frame(
  Sample = c('Control_1', 'Control_2', 'Treatment_1', 'Treatment_2'),
  Group  = c('Control', 'Control', 'Treatment', 'Treatment')
)

# Step 1: Normalize data
normalized_obj <- proteo.normalize(raw_data, metadata)

# Step 2: Run enrichment (min_genes = 1 for small example)
proteo.enrich(normalized_obj, min_genes = 1)
")
    return(invisible(NULL))
  }

  db <- match.arg(db)
  orgdb_pkg <- switch(db,
                      hs = "org.Hs.eg.db",
                      mm = "org.Mm.eg.db")

  # ---- Check required packages ----
  pkgs <- c("clusterProfiler", orgdb_pkg, "ggplot2")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse=", "))

  # ---- Load OrgDb ----
  OrgDb <- getExportedValue(orgdb_pkg, orgdb_pkg)

  # ---- Helper: clean UniProt IDs ----
  .clean_uniprot <- function(x) {
    x <- sub("^\\w+\\|", "", x)
    x <- sub("\\|.*$", "", x)
    x <- sub("-\\d+$", "", x)
    trimws(x)
  }

  # ---- Extract presence_per_group ----
  if(is.list(normalized_data) && !is.null(normalized_data$presence_per_group)) {
    presence_per_group <- normalized_data$presence_per_group
  } else if(!is.null(attr(normalized_data, "presence_per_group"))) {
    presence_per_group <- attr(normalized_data, "presence_per_group")
  } else if(is.data.frame(normalized_data) && "ProteinID" %in% colnames(normalized_data)) {
    presence_per_group <- data.frame(ProteinID = normalized_data$ProteinID, ALL = 1, check.names = FALSE)
  } else {
    stop("The object does not contain 'presence_per_group' nor a 'ProteinID' column.")
  }

  presence_per_group <- as.data.frame(presence_per_group)
  if(!"ProteinID" %in% colnames(presence_per_group)) stop("'presence_per_group' must contain a 'ProteinID' column.")

  # ---- Groups and proteins ----
  groups <- setdiff(colnames(presence_per_group), "ProteinID")
  if(length(groups) == 0) stop("No groups found.")

  protein_ids <- .clean_uniprot(presence_per_group$ProteinID)
  proteins_by_group <- lapply(groups, function(g) protein_ids[presence_per_group[[g]] > 0])
  names(proteins_by_group) <- groups
  proteins_unique <- lapply(groups, function(g) setdiff(proteins_by_group[[g]], unlist(proteins_by_group[groups != g])))
  names(proteins_unique) <- groups
  proteins_shared <- Reduce(intersect, proteins_by_group)

  # ---- Messages ----
  if (verbose) {
    report <- paste0(
      "\n=== Unique proteins per group ===\n",
      paste0(
        vapply(groups,
               function(grp) sprintf("%s: %d unique proteins\n", grp, length(proteins_unique[[grp]])),
               FUN.VALUE = character(1)),
        collapse = ""
      ),
      sprintf("Proteins shared across all groups: %d\n", length(proteins_shared)),
      if (!results) "\nTo see detailed results, use 'results = TRUE'.\n" else ""
    )
    message(report)
  }

  # ---- GO enrichment ----
  ego_list <- list()
  proteins_used <- list()

  for(grp in groups) {

    # Select proteins: unique or all
    genes <- if(proteins == "unique") proteins_unique[[grp]] else proteins_by_group[[grp]]
    proteins_used[[grp]] <- genes

    if(length(genes) < min_genes) {
      message("\nInsufficient number of proteins in ", grp, " (minimum: ", min_genes, ")")
      next
    }

    # Map UniProt -> Entrez
    entrez_ids <- clusterProfiler::bitr(genes, fromType="UNIPROT", toType="ENTREZID", OrgDb=OrgDb)
    if(nrow(entrez_ids) == 0) {
      message("\nNo genes mapped to EntrezID for ", grp, ". Cleaned IDs: ", paste(head(genes,5), collapse=", "))
      next
    }

    # GO enrichment
    ego <- clusterProfiler::enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = OrgDb,
      keyType = "ENTREZID",
      ont = ont,
      readable = TRUE,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )

    # Safe check: enrichResult or data.frame
    ego_res <- if(inherits(ego, "enrichResult")) ego@result else as.data.frame(ego)

    # Filter by FDR
    ego_res <- ego_res[ego_res$p.adjust <= p_adjust_cutoff, ]

    if(nrow(ego_res) > 0 && inherits(ego, "enrichResult")) {
      ego_list[[grp]] <- ego
      if(show_plot) {
        message("\n=== Dotplot GO ", ont, " - group ", grp, " ===")
        print(clusterProfiler::dotplot(ego, showCategory=15) +
                ggplot2::ggtitle(paste0("GO ", ont, " - ", grp)))
      }
    } else {
      message("\nNo significant GO terms (p.adjust < ", p_adjust_cutoff, ") for ", grp)
    }
  }

  if(results) return(list(ego_list=ego_list, proteins_used=proteins_used))
  invisible(NULL)
}
