#' Print Method for proteoNorm Objects
#'
#' Custom print method for objects of class \code{proteoNorm}.
#'
#' @param x An object of class \code{proteoNorm}.
#' @param ... Additional arguments (ignored).
#'
#' @export

print.proteoNorm <- function(x, ...) {
  message("================================")
  message("Object of class 'proteoNorm'")
  message("================================")

  df <- x$normalized_matrix
  if (is.null(df)) {
    message("No normalized matrix found in the object.")
    message("================================")
    return(invisible(x))
  }

  # Identifica coluna de IDs automaticamente (não numérica)
  id_col <- which(!sapply(df, is.numeric))

  # Conta amostras corretamente
  sample_count <- ncol(df) - length(id_col)

  message("Samples:  ", sample_count)
  message("Proteins: ", nrow(df))
  message("Method:   ", x$method)

  if (!is.null(x$checks) && length(x$checks) > 0) {
    message("--- Checks ---")
    for (c in x$checks) message(c)
  }

  message("================================")
  invisible(x)
}




