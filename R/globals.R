# R/globals.R
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "ProteinID", "Sample", "Group", "Abundance", "htest", "p_value", "signif_label",
    "value", "variable", "n", "x", "y", "r_val", "p_val"
  ))
}
