#' Normalize Proteomics Data
#'
#' Performs normalization on raw proteomics abundance data. Supports log transformation,
#' optional outlier removal, global imputation for missing values, and filtering
#' proteins with excessive missing values per group.
#'
#' @param raw_data A \code{data.frame} with protein abundances. Must contain a column with
#'   protein IDs (e.g., "ProteinID") and columns corresponding to sample names.
#' @param metadata A \code{data.frame} with sample metadata. Must contain columns for
#'   sample names and group/condition labels.
#' @param method Character. Normalization method. Options: \code{"log2"} (default), \code{"log10"}, or \code{"none"}.
#' @param remove_zeros Logical. Whether to remove proteins with zero values in all samples. Default: \code{FALSE}.
#' @param remove_outliers Logical. Whether to replace extreme values with \code{NA}. Default: \code{TRUE}.
#' @param impute Logical. Whether to perform global imputation for missing/zero values. Default: \code{TRUE}.
#' @param impute_shift Numeric. Shift factor for global imputation mean. Default: 1.8.
#' @param impute_scale Numeric. Scale factor for global imputation standard deviation. Default: 0.3.
#' @param max_missing_prop_per_group Numeric. Maximum proportion of missing/zero values allowed per group. Default: 0.5.
#' @param pseudo_count Numeric. Small value added before log transformation. Default: 1e-6.
#' @param clip_negative_log Logical. Whether to set negative values to 0 after log transformation. Default: \code{TRUE}.
#' @param assign_result Logical. Whether to assign the resulting \code{proteoNorm} object in the calling environment. Default: \code{TRUE}.
#' @param assign_name Character. Name to assign the resulting object. Default: \code{"normalized_data"}.
#' @param envir Environment in which to assign the result if \code{assign_result = TRUE}. Default: parent.frame().
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return An object of class \code{proteoNorm}, which is a list containing:
#'   \itemize{
#'     \item \code{<method>_normalized}: Normalized data.frame with \code{ProteinID} and samples.
#'     \item \code{raw_preserved}: Original raw data for reference.
#'     \item \code{presence_per_group}: Matrix indicating protein presence/absence per group.
#'     \item \code{method}: Normalization method used.
#'     \item \code{checks}: List of checks and actions applied during normalization.
#'   }
#'
#' @examples
#' raw_data <- data.frame(
#'   ProteinID = c("P1", "P2", "P3", "P4"),
#'   Control_1 = c(1.2e6, 0, 5e5, 2e5),
#'   Control_2 = c(1.1e6, 0, 4.8e5, 2.2e5),
#'   Treatment_1 = c(0.9e6, 0, 5.2e5, 1.8e5),
#'   Treatment_2 = c(1.0e6, 0, 5.1e5, 2.0e5)
#' )
#' metadata <- data.frame(
#'   Sample = c("Control_1", "Control_2", "Treatment_1", "Treatment_2"),
#'   Group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' normalized_data <- proteo.normalize(raw_data, metadata)
#' str(normalized_data)
#'
#' @export
#'
#' @references
#' Pearson, K. 1901. On lines and planes of closest fit to systems of points in space.
#'   \emph{Philosophical Magazine} 2(11): 559–572. <doi:10.1080/14786440109462720>

# --- Main function ---
proteo.normalize <- function(raw_data, metadata, method = "log2",
                             remove_zeros = FALSE, remove_outliers = TRUE,
                             impute = TRUE, impute_shift = 1.8, impute_scale = 0.3,
                             max_missing_prop_per_group = 0.5, pseudo_count = 1e-6,
                             clip_negative_log = TRUE,
                             assign_result = TRUE, assign_name = "normalized_data",
                             envir = parent.frame(),
                             help = FALSE) {

  # --- Help message ---
  if (help || missing(raw_data) || missing(metadata)) {
    message("
Function proteo.normalize()

Description:
  Performs normalization of raw proteomics abundance data.
  Steps include optional log transformation, outlier removal, global imputation
  for missing values, and filtering proteins with excessive missing values per group.
  Returns an object of class 'proteoNorm'.

Usage:
  proteo.normalize(raw_data, metadata, method = 'log2', ...)

Main arguments:
  raw_data                   Data frame with protein abundances. Must include a protein ID column.
  metadata                   Data frame with sample metadata. Must include sample names and group labels.
  method                     Normalization method: 'log2' (default), 'log10', or 'none'.
  remove_zeros               Logical. Remove proteins with zero values in all samples? Default FALSE.
  remove_outliers             Logical. Replace extreme values with NA? Default TRUE.
  impute                     Logical. Perform global imputation for missing/zero values? Default TRUE.
  impute_shift               Numeric. Shift factor for imputation mean. Default 1.8.
  impute_scale               Numeric. Scale factor for imputation SD. Default 0.3.
  max_missing_prop_per_group Numeric. Maximum allowed missing proportion per group. Default 0.5.
  pseudo_count               Numeric. Small value added before log transformation. Default 1e-6.
  clip_negative_log          Logical. Set negative values to 0 after log? Default TRUE.
  assign_result              Logical. Assign result in calling environment? Default TRUE.
  assign_name                Character. Name of assigned object. Default 'normalized_data'.
  envir                      Environment for assignment if assign_result = TRUE.
  help                       Logical. If TRUE, prints this help message.

Return:
  An object of class 'proteoNorm' containing:
    - <method>_normalized : Normalized data frame with ProteinID and samples
    - raw_preserved       : Original raw data for reference
    - presence_per_group  : Matrix indicating protein presence/absence per group
    - method              : Normalization method used
    - checks              : List of validation steps and actions performed

Example:
  raw_data <- data.frame(
    ProteinID = c('P1','P2','P3','P4'),
    Control_1 = c(1.2e6,0,5e5,2e5),
    Control_2 = c(1.1e6,0,4.8e5,2.2e5),
    Treatment_1 = c(0.9e6,0,5.2e5,1.8e5),
    Treatment_2 = c(1.0e6,0,5.1e5,2.0e5)
  )

  metadata <- data.frame(
    Sample = c('Control_1','Control_2','Treatment_1','Treatment_2'),
    Group  = c('Control','Control','Treatment','Treatment')
  )

  normalized_data <- proteo.normalize(raw_data, metadata)
  print(normalized_data)
")
    return(invisible(NULL))
  }

  checks <- list()

  # --- Validations ---
  if (!is.data.frame(raw_data)) stop("'raw_data' must be a data.frame.")
  if (!is.data.frame(metadata)) stop("'metadata' must be a data.frame.")

  # ProteinID
  col_id <- which(names(raw_data) %in% c("ProteinID","Protein_ID","Protein",
                                         "Accession","Protein.Accession","protein_id"))
  if (length(col_id) != 1) stop("Could not detect Protein ID column.")
  checks$protein_id <- paste("Protein ID column detected:", names(raw_data)[col_id])

  # Sample and Group
  col_sample <- which(names(metadata) %in% c("Sample","sample","SAMPLE"))
  col_group <- which(names(metadata) %in% c("Group","group","GROUP",
                                            "Condition","condition","COND"))
  if (length(col_sample) != 1) stop("Sample column not detected in metadata.")
  if (length(col_group) != 1) stop("Group column not detected in metadata.")

  samples <- metadata[[col_sample]]
  groups <- metadata[[col_group]]

  expr_data <- raw_data[, samples, drop = FALSE]

  # --- Remove zero proteins ---
  if (remove_zeros) {
    nonzero_rows <- rowSums(expr_data != 0, na.rm = TRUE) > 0
    n_removed <- sum(!nonzero_rows)
    expr_data <- expr_data[nonzero_rows, , drop = FALSE]
    raw_data <- raw_data[nonzero_rows, , drop = FALSE]
    checks$removed_all_zeros <- paste("Proteins removed (all zero):", n_removed)
  }

  # --- Remove outliers ---
  if (remove_outliers) {
    remove_outlier_row <- function(x) {
      qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- qnt[2] - qnt[1]
      x[x < qnt[1] - 3*iqr | x > qnt[2] + 3*iqr] <- NA
      x
    }
    expr_data <- t(apply(expr_data, 1, remove_outlier_row))
    checks$outliers <- "Outliers replaced with NA"
  }

  # --- Presence/absence per group (antes do filtro de missing) ---
  if (nrow(expr_data) > 0) {
    presence_per_group <- t(apply(expr_data, 1, function(x)
      tapply(x, groups, function(vals) as.integer(any(vals > 0)))
    ))
    presence_per_group <- cbind(ProteinID = raw_data$ProteinID, presence_per_group)
  } else {
    presence_per_group <- data.frame()
  }

  # --- Filter proteins with too many missing values ---
  if (nrow(expr_data) > 0) {
    remove_idx <- sapply(1:nrow(expr_data), function(i) {
      by_group <- tapply(expr_data[i, ], groups, function(vals) mean(is.na(vals) | vals == 0))
      all(by_group > max_missing_prop_per_group)
    })
    expr_data <- expr_data[!remove_idx, , drop = FALSE]
    raw_data <- raw_data[!remove_idx, , drop = FALSE]
    if (nrow(presence_per_group) > 0) presence_per_group <- presence_per_group[!remove_idx, , drop = FALSE]
    checks$missing_filter <- paste("Proteins removed due to excessive missing values:", sum(remove_idx))
  }

  # --- Stop if no proteins remain ---
  if (nrow(expr_data) == 0) {
    warning("No proteins remain after filtering. Returning empty normalized object.")
    norm_matrix <- matrix(nrow = 0, ncol = length(samples))
    colnames(norm_matrix) <- samples
  } else {
    # --- Normalization ---
    norm_matrix <- as.matrix(expr_data)
    norm_matrix[!is.finite(norm_matrix) | norm_matrix < 0] <- 0
    if (method == "log2") norm_matrix <- log2(norm_matrix + pseudo_count)
    if (method == "log10") norm_matrix <- log10(norm_matrix + pseudo_count)

    # --- Global imputation ---
    if (impute) {
      observed <- norm_matrix[norm_matrix > 0 & is.finite(norm_matrix)]
      imp_mean <- median(observed, na.rm = TRUE) - impute_shift * sd(observed, na.rm = TRUE)
      imp_sd <- sd(observed, na.rm = TRUE) * impute_scale
      miss_idx <- which(norm_matrix <= 0 | is.na(norm_matrix))
      if (length(miss_idx) > 0) norm_matrix[miss_idx] <- rnorm(length(miss_idx), mean = imp_mean, sd = imp_sd)
      checks$imputation <- sprintf("Global imputation applied (mean = %.2f, sd = %.2f)", imp_mean, imp_sd)
    }

    # --- Clip negatives ---
    if (clip_negative_log) {
      n_neg <- sum(norm_matrix < 0, na.rm = TRUE)
      if (n_neg > 0) norm_matrix[norm_matrix < 0] <- 0
      checks$clip_neg <- paste("Negative values after log adjusted to 0:", n_neg)
    }
  }

  # --- Create final data.frames ---
  norm_name <- paste0(method, "_normalized")
  normalized_df <- cbind(ProteinID = raw_data$ProteinID, as.data.frame(norm_matrix))
  raw_preserved <- cbind(ProteinID = raw_data$ProteinID, raw_data[, samples, drop = FALSE])

  # --- Build final object ---
  obj <- list()
  obj[[norm_name]] <- normalized_df
  obj$raw_preserved <- raw_preserved
  obj$presence_per_group <- presence_per_group
  obj$method <- method
  obj$checks <- checks
  class(obj) <- "proteoNorm"

  # --- Assign and print ---
  if (assign_result) {
    assign(assign_name, obj, envir = envir)
    print(obj)
    invisible(obj)
  } else {
    return(obj)
  }
}
