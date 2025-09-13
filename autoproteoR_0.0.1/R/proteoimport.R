#' Import proteomics data
#'
#' Validates and prepares raw proteomics data for downstream analysis. Detects
#' the protein ID column, checks quantification columns, and optionally validates
#' metadata. Returns an object of class `proteoR`.
#'
#' @param raw_data Data frame or tibble containing raw protein abundances. Must
#'   have a protein identifier column.
#' @param metadata Optional data frame containing sample metadata. Must contain
#'   a column named `Sample` if provided.
#' @param strict Logical; if TRUE (default), the function stops on validation
#'   errors. If FALSE, warnings are issued instead.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return An object of class `proteoR` containing:
#'   \describe{
#'     \item{data}{The raw data as provided.}
#'     \item{metadata}{Metadata if provided; otherwise NULL.}
#'     \item{protein_id_col}{Index of the column identified as protein ID.}
#'     \item{checks}{List of validation checks performed (including warnings).}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # --- Small example ---
#' raw_data <- data.frame(
#'   ProteinID   = c("P001", "P002"),
#'   Control_1   = c(10000, 20000),
#'   Treatment_1 = c(8000, 19000)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control_1", "Treatment_1"),
#'   Group  = c("Control", "Treatment")
#' )
#'
#' obj <- proteo.import(raw_data, metadata)
#' print(obj)
#' }


proteo.import <- function(raw_data, metadata = NULL, strict = TRUE, help = FALSE) {

  # --- Help message ---
  if (help || missing(raw_data)) {
    message("
Function proteo.import()

Description:
  Imports and validates raw proteomics data for downstream analysis.
  Detects the protein identifier column, checks quantification columns,
  and optionally validates sample metadata. Returns an object of class 'proteoR'.

Usage:
  proteo.import(raw_data, metadata = NULL, strict = TRUE)

Arguments:
  raw_data  Data frame or tibble with raw protein abundances. Must include a protein ID column.
  metadata  Optional data frame with sample metadata. Must have a 'Sample' column if provided.
  strict    Logical. If TRUE (default), stops on validation errors. If FALSE, issues warnings.
  help      Logical. If TRUE, prints this help message.

Return:
  An object of class 'proteoR' containing:
    - data          The raw data as provided.
    - metadata      Metadata if provided; otherwise NULL.
    - protein_id_col Index of the protein ID column detected.
    - checks        List of validation checks performed (including warnings).

Example:
  # Small example
  raw_data <- data.frame(
    ProteinID   = c('P001', 'P002'),
    Control_1   = c(10000, 20000),
    Treatment_1 = c(8000, 19000)
  )

  metadata <- data.frame(
    Sample = c('Control_1', 'Treatment_1'),
    Group  = c('Control', 'Treatment')
  )

  obj <- proteo.import(raw_data, metadata)
  print(obj)
")
    return(invisible(NULL))
  }
  checks <- list()
  errors <- character()
  warnings <- character()

  add_error <- function(text) {
    errors <<- c(errors, text)
  }
  add_warning <- function(text) {
    warnings <<- c(warnings, text)
  }

  ## 1. Object type
  if (!is.data.frame(raw_data)) {
    add_error("'raw_data' must be a data.frame or tibble.")
  } else {
    checks$raw_data <- "Object type: OK"
  }

  ## 2. Detect Protein ID column
  col_id <- 1 # default
  possible_names <- c("ProteinID", "Protein_ID", "Protein", "Accession", "Protein.Accession", "protein_id")

  match_idx <- match(tolower(possible_names), tolower(colnames(raw_data)))
  match_idx <- match_idx[!is.na(match_idx)]

  if (length(match_idx) > 0) {
    col_id <- match_idx[1]
    checks$protein_id <- paste0("Protein ID column detected: '", colnames(raw_data)[col_id], "'")
  } else {
    add_warning("Protein ID column not detected; assuming the 1st column.")
    checks$protein_id <- "Protein ID column not detected; assuming the 1st column."
  }

  # Must be text
  if (!is.character(raw_data[[col_id]])) {
    add_error("The Protein ID column does not appear to be text. Please check your data.")
  }

  ## 3. Quantification columns
  quant_cols <- setdiff(seq_len(ncol(raw_data)), col_id)

  if (!all(sapply(raw_data[quant_cols], is.numeric))) {
    add_error("All quantification columns must be numeric.")
  } else {
    checks$quant_cols <- "Quantification columns: OK"
  }

  ## 4. Optional metadata
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      add_error("'metadata' must be a data.frame or tibble.")
    } else {
      if (!"Sample" %in% colnames(metadata)) {
        add_error("Metadata must contain a column named 'Sample'.")
      } else {
        if (!all(metadata$Sample %in% colnames(raw_data)[quant_cols])) {
          add_warning("Some samples in metadata are not present in 'raw_data'.")
        }
      }
    }
    checks$metadata <- "Metadata provided and validated"
  } else {
    checks$metadata <- "No metadata provided"
  }

  # Stop if errors and strict = TRUE
  if (length(errors) > 0) {
    if (strict) {
      stop(paste("Data validation error(s):\n-", paste(errors, collapse = "\n- ")), call. = FALSE)
    } else {
      warnings <- c(warnings, errors)
    }
  }

  # Only warnings, no errors
  if (length(errors) == 0 && length(warnings) > 0) {
    obj <- list(
      data = NULL,
      metadata = NULL,
      protein_id_col = NULL,
      checks = list(warning = paste(warnings, collapse = "; "))
    )
    class(obj) <- "proteoR"
    return(obj)
  }

  # No errors/warnings → return full object
  obj <- list(
    data = raw_data,
    metadata = metadata,
    protein_id_col = col_id,
    checks = checks
  )

  class(obj) <- "proteoR"
  return(obj)
}

# --- S3 print method ---
#' @exportS3Method print proteoR
print.proteoR <- function(x, ...) {
  cat("Object of class 'proteoR'\n")
  cat("------------------------------\n")
  if (!is.null(x$checks$warning)) {
    cat("Warning:", x$checks$warning, "\n")
  } else {
    for (check in x$checks) {
      cat(check, "\n")
    }
    cat("------------------------------\n")
    cat("Proteins:", ifelse(is.null(x$data), "NA", nrow(x$data)), "\n")
    cat("Samples:", ifelse(is.null(x$data), "NA", ncol(x$data) - 1), "\n")
    if (!is.null(x$metadata)) {
      cat("Metadata with", nrow(x$metadata), "rows\n")
    }
  }
  invisible(x)
}











