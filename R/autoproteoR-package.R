#' autoproteoR: Quick Reference and Main Functions for Proteomics Analysis
#'
#' autoproteoR provides a complete suite of functions for proteomics data analysis,
#' from raw data import to exploratory and statistical analysis. This help topic
#' serves as a quick reference and guide for users.
#'
#' Main Functions Overview
#' | Function | Description |
#' |----------|-------------|
#' | `proteo.import()`  | Import and validate raw proteomics data |
#' | `proteo.normalize()`| Normalize and filter data |
#' | `proteo.qc()`      | QC plots (boxplot, PCA) |
#' | `proteo.corr()`    | Compute sample correlations |
#' | `proteo.dimred()`  | PCA and UMAP dimensionality reduction |
#' | `proteo.enrich()`  | GO enrichment analysis |
#' | `proteo.heatmap()` | For top variable proteins |
#' | `proteo.ttest()`   | T-test or Mann-Whitney test |
#' | `proteo.venn()`    | Venn diagrams |
#' | `proteo.volcano()` | For differentially expressed proteins |
#' | `proteo.roc()`     | ROC curve for classification analysis |
#'
#'
#' Contact and Contributions
#' For suggestions, bug reports, or contributions, see the
#' [GitHub repository](https://github.com/Luiz-Garcia-R/autoproteoR).
#'
#'
#' Example workflow for proteo.import and proteo.normalize
#'
#' A small synthetic proteomics dataset with 50 proteins and 5 samples per group.
#' This dataset is for demonstration purposes only.
#'
#' @examples
#'
#' proteins <- paste0("P", sprintf("%03d", 1:50))
#'
#' raw_data <- data.frame(
#'   ProteinID = proteins,
#'   Control_1   = rnorm(50, mean = 10000, sd = 2000),
#'   Control_2   = rnorm(50, mean = 10000, sd = 2000),
#'   Control_3   = rnorm(50, mean = 10000, sd = 2000),
#'   Control_4   = rnorm(50, mean = 10000, sd = 2000),
#'   Control_5   = rnorm(50, mean = 10000, sd = 2000),
#'   Treatment_1 = rnorm(50, mean = 20000, sd = 2500),
#'   Treatment_2 = rnorm(50, mean = 20000, sd = 2500),
#'   Treatment_3 = rnorm(50, mean = 20000, sd = 2500),
#'   Treatment_4 = rnorm(50, mean = 30000, sd = 2500),
#'   Treatment_5 = rnorm(50, mean = 30000, sd = 2500)
#' )
#'
#' metadata <- data.frame(
#'   Sample = colnames(raw_data)[-1],
#'   Group  = rep(c("Control", "Treatment"), each = 5)
#' )
#'
#'
#' obj <- proteo.import(raw_data, metadata)
#' normalized <- proteo.normalize(raw_data, metadata)
#'
#' proteo.qc(normalized, metadata)
#'
#' corr_mat <- proteo.corr(normalized, metadata)
#' dimred_res <- proteo.dimred(normalized, metadata)
#'
#' proteo.volcano(normalized, metadata,
#'                group_col = "Group",
#'                padj_threshold = 0.05,
#'                log2fc_threshold = 1,
#'                results = TRUE)
#'
#' @name autoproteoR
#'
"_PACKAGE"


