# R/imports.R

# --- Funções importadas ---
#' @importFrom stats sd aggregate shapiro.test ks.test t.test chisq.test var.test setNames cor cor.test median quantile rnorm
#' @importFrom utils globalVariables head
#' @importFrom dplyr %>% left_join group_by summarise select everything
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_point geom_smooth annotate theme_minimal
#'   theme element_text labs scale_y_continuous
#' @importFrom pheatmap pheatmap
#' @importFrom ggvenn ggvenn
#' @importFrom RColorBrewer brewer.pal

# --- Variáveis globais ---
if(getRversion() >= "2.15.1")  utils::globalVariables(
  c("PC1", "PC2", "UMAP1", "UMAP2", "Expression", "log2FoldChange",
    "pvalue", "group_color", "Protein")
)

NULL
