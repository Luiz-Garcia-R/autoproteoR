#' Compute and visualize correlation of proteomics data
#'
#' Computes pairwise correlations among samples or groups for proteomics data.
#' Automatically chooses Pearson or Spearman correlation based on normality and ties,
#' and provides heatmap and scatterplot visualizations.
#'
#' @param normalized_data A \code{data.frame} or \code{proteonorm} object containing normalized protein expression.
#' @param group_level Logical. If TRUE (default), correlation is computed at the group level (averaging samples per group).
#' @param method Character. Correlation method: "pearson", "spearman", "kendall", or "auto" (default).
#' @param plot Logical. If TRUE (default), a heatmap is plotted, and a scatterplot is drawn if only 2 columns.
#' @param help Logical. If TRUE, prints this help message and example usage. Default = FALSE.
#'
#' @return A correlation matrix.
#'
#' @examples
#' # Small example: raw data
#' raw_data <- data.frame(
#'   ProteinID = paste0("P", 1:5),
#'   Control1 = c(100, 200, 150, 300, 250),
#'   Control2 = c(110, 210, 160, 310, 260),
#'   Treatment1 = c(300, 100, 200, 150, 250),
#'   Treatment2 = c(310, 90, 210, 140, 260)
#' )
#'
#' metadata <- data.frame(
#'   Sample = c("Control1", "Control2", "Treatment1", "Treatment2"),
#'   Group = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' # Step 1: Normalize raw data
#' normalized_obj <- proteo.normalize(raw_data, metadata)
#'
#' # Step 2: Compute correlation matrix
#' corr_mat <- proteo.corr(normalized_obj, metadata)
#'
#' @references
#' Pearson correlation: Pearson 1895 <doi:10.1098/rspl.1895.0041>
#' Spearman correlation: Spearman 1904 <doi:10.1037/1082-989X.9.2.220>
#'
#' @export


proteo.corr <- function(normalized_data, group_level = TRUE,
                        method = "auto", plot = TRUE, help = FALSE) {

  if (help || missing(normalized_data)) {
    message("
Function proteo.corr()

Description:
  Computes pairwise correlations among samples or groups for normalized proteomics data.
  Automatically selects Pearson, Spearman or Kendall if 'auto', and provides heatmap
  and a styled scatterplot when there are exactly 2 columns to compare.

Usage:
  proteo.corr(normalized_data, group_level = TRUE, method = 'auto', plot = TRUE)
")
    return(invisible(NULL))
  }

  # --- Extract metadata and matrix if proteonorm ---
  if ("proteonorm" %in% class(normalized_data)) {
    metadata <- normalized_data$metadata
    expr_mat <- as.data.frame(normalized_data$normalized_matrix, stringsAsFactors = FALSE)
  } else {
    expr_mat <- normalized_data
    metadata <- NULL
  }

  # --- Detect ProteinID column and prepare expression matrix ---
  possible_id_names <- c("ProteinID", "Protein_ID", "Protein", "Accession", "Protein.Accession", "protein_id")
  col_id <- which(names(expr_mat) %in% possible_id_names)
  if (length(col_id) == 1) {
    protein_names <- as.character(expr_mat[[col_id]])
    protein_names[is.na(protein_names) | protein_names == ""] <- "Protein_Unknown"
    protein_names <- make.unique(protein_names)
    expr_mat <- expr_mat[, -col_id, drop = FALSE]
    rownames(expr_mat) <- protein_names
  } else {
    # if no ProteinID column, try to use rownames
    if (is.null(rownames(expr_mat)) || any(rownames(expr_mat) == "")) {
      stop("Could not detect ProteinID column and rownames are not set.")
    }
  }

  # Convert to data.frame (safe) and ensure numeric columns
  expr_df <- as.data.frame(expr_mat, stringsAsFactors = FALSE)
  # Ensure numeric columns (coerce if necessary)
  expr_df[] <- lapply(expr_df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (all(is.na(as.numeric(as.character(col))))) return(as.numeric(rep(NA, length(col))))
    as.numeric(as.character(col))
  })

  # --- Aggregate by group if requested ---
  if (group_level) {
    if (is.null(metadata)) stop("metadata must be provided (or pass a proteonorm object) when group_level = TRUE")
    if (!all(c("Sample","Group") %in% colnames(metadata))) stop("Metadata must contain 'Sample' and 'Group' columns")

    df_long <- tidyr::pivot_longer(
      dplyr::mutate(expr_df, ProteinID = rownames(expr_df)),
      -ProteinID, names_to = "Sample", values_to = "Abundance"
    )
    df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

    df_grouped <- dplyr::summarise(
      dplyr::group_by(df_long, ProteinID, Group),
      Abundance = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

    df_grouped <- tidyr::pivot_wider(df_grouped, names_from = Group, values_from = Abundance)
    corr_df <- as.data.frame(df_grouped[, -1, drop = FALSE])
  } else {
    corr_df <- expr_df
  }

  # Remove columns that are all NA or non-numeric (safety)
  ok_cols <- vapply(corr_df, function(col) is.numeric(col) && any(!is.na(col)), logical(1))
  if (!all(ok_cols)) corr_df <- corr_df[, ok_cols, drop = FALSE]
  if (ncol(corr_df) < 1) stop("No numeric columns available for correlation after preprocessing.")

  # --- Choose correlation method (auto logic kept simple here) ---
  if (method == "auto") {
    # Quick normality probe per column using Shapiro when small; fallback to AD/KS if larger
    aplicar_teste_normalidade <- function(x) {
      x <- x[is.finite(x)]
      n <- length(x)
      if (n < 3) return(list(p = NA, metodo = "none"))
      if (n <= 50) {
        p <- tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA)
        metodo <- "Shapiro-Wilk"
      } else if (n <= 300) {
        if (!requireNamespace("nortest", quietly = TRUE)) {
          p <- NA; metodo <- "no-nortest"
        } else {
          p <- tryCatch(nortest::ad.test(x)$p.value, error = function(e) NA)
          metodo <- "Anderson-Darling"
        }
      } else {
        p <- tryCatch(stats::ks.test(x, "pnorm", mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))$p.value, error = function(e) NA)
        metodo <- "Kolmogorov-Smirnov"
      }
      list(p = p, metodo = metodo)
    }

    normal_res <- lapply(as.data.frame(corr_df), aplicar_teste_normalidade)
    normality <- vapply(normal_res, function(z) if (is.na(z$p)) FALSE else z$p > 0.05, logical(1))
    # detect ties
    has_ties <- any(vapply(corr_df, function(x) anyDuplicated(x) > 0, logical(1)))
    # detect outliers proportion (IQR rule) across both columns when there are 2 columns
    detect_outliers_prop <- function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 3) return(0)
      q1 <- stats::quantile(v, 0.25, na.rm = TRUE); q3 <- stats::quantile(v, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      if (iqr == 0) return(0)
      mean(v < (q1 - 1.5*iqr) | v > (q3 + 1.5*iqr), na.rm = TRUE)
    }
    out_prop <- if (ncol(corr_df) >= 2) max(vapply(corr_df, detect_outliers_prop, numeric(1))) else 0

    metodo_usado <- if (out_prop > 0.05) {
      "kendall"
    } else if (all(normality) && !has_ties) {
      "pearson"
    } else {
      "spearman"
    }

    message("Automatically selected correlation method: ", metodo_usado)
  } else {
    metodo_usado <- method
  }

  # --- Correlation matrix ---
  corr_mat <- stats::cor(corr_df, method = metodo_usado, use = "pairwise.complete.obs")

  # --- Create "r = ..." matrix for display ---
  corr_text <- matrix(paste0("r = ", round(corr_mat, 2)),
                      nrow = nrow(corr_mat), ncol = ncol(corr_mat),
                      dimnames = dimnames(corr_mat))

  # --- Plot heatmap + scatterplot (improved) ---
  if (plot) {
    pheatmap::pheatmap(
      corr_mat, display_numbers = corr_text, number_color = "black",
      main = paste("Correlation (", metodo_usado, ")"),
      cluster_rows = FALSE, cluster_cols = FALSE
    )

    # If exactly 2 columns, produce styled scatterplot (estilo 1)
    if (ncol(corr_df) == 2) {
      cols <- colnames(corr_df)
      # prepare data
      dados <- data.frame(x = corr_df[[1]], y = corr_df[[2]])
      nome_x <- cols[1]; nome_y <- cols[2]

      # perform cor.test using chosen method
      teste <- tryCatch(stats::cor.test(dados$x, dados$y, method = metodo_usado), error = function(e) NULL)
      if (!is.null(teste)) {
        coef_val <- round(unname(teste$estimate), 3)
        pval <- teste$p.value
        p_texto <- if (is.na(pval)) "p = NA" else if (pval < 0.001) "p < 0.001" else paste0("p = ", signif(pval, 3))
      } else {
        coef_val <- NA; p_texto <- "p = NA"
      }

      interprete <- if (!is.na(coef_val)) {
        cut(abs(coef_val),
            breaks = c(-Inf, 0.3, 0.5, 0.7, 0.9, Inf),
            labels = c("Very weak", "Weak", "Moderated", "Strong", "Very strong"),
            right = FALSE)
      } else {
        NA
      }

      simbolo <- ifelse(tolower(metodo_usado) == "kendall", "tau", "r")
      subtitulo <- paste0(simbolo, " = ", coef_val, " | ", p_texto, " | ", interprete)

      # choose smoothing method
      smooth_method <- if (tolower(metodo_usado) == "pearson") "lm" else "loess"

      # estilo 1 aesthetics
      g <- ggplot2::ggplot(dados, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(alpha = 0.60, size = 2.5, col = "grey30") +
        ggplot2::geom_smooth(method = smooth_method, se = FALSE, linetype = "dashed", color = "steelblue") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(title = paste(tools::toTitleCase(metodo_usado), "Correlation"),
                      subtitle = subtitulo, x = nome_x, y = nome_y)

      print(g)
    }
  }

  invisible(corr_mat)
}
