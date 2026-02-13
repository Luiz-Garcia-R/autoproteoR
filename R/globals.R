if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(

    # Identificadores / dados
    "ProteinID", "Protein", "Sample", "Group", "Abundance",
    "Expression", "value", "variable",

    # Estatística
    "p_value", "p_val", "pvalue", "r_val", "n", "x", "y", "htest",

    # Fold-change / ranking
    "logFC", "log2FoldChange", "direction",
    "Metric", "Rank",

    # Plot / visual
    "PlotLabel", "signif_label", "group_color",

    # Dimensionalidade
    "PC1", "PC2", "UMAP1", "UMAP2"

  ))
}
