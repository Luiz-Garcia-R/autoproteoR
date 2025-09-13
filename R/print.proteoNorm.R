#' Print Method for proteoNorm Objects
#'
#' Custom print method for objects of class \code{proteoNorm}.
#'
#' @param x An object of class \code{proteoNorm}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.proteoNorm <- function(x, ...) {
  cat("Object of class 'proteoNorm'\n")
  cat("------------------------------\n")
  for (check in x$checks) cat(check, "\n")
  cat("------------------------------\n")
  cat("Normalization method:", x$method, "\n")
  norm_df <- x[[grep("_normalized$", names(x))]]
  cat("Proteins:", nrow(norm_df), "\n")
  cat("Samples:", ncol(norm_df) - 1, "\n")
  invisible(x)
}
