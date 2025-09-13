# autoproteoR

<!-- badges: start -->
<!-- badges: end -->

**autoproteoR** is an R package designed to streamline proteomics data analysis, from raw data import to exploratory and statistical evaluation.  
It provides functions for normalization, quality control, correlation, dimensionality reduction, enrichment, heatmaps, volcano plots, statistical testing, and Venn diagrams.

**Note:** The package is primarily optimized for **pairwise analyses** (e.g., Control vs. Treatment).  
Workflows involving multiple groups can be explored, but most functions are tuned for two-group comparisons.

## Installation

You can install the development version directly from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install autoproteoR from GitHub
devtools::install_github("https://github.com/Luiz-Garcia-R/autoproteoR.git")
```

## Input Data Format (Very Important)

To use autoproteoR, you need two input data frames:
  1. raw_data – protein quantifications (output from MaxQuant, PatternLab, FragPipe, etc.)
    - Must contain a first column called ProteinID
    - Other columns must be samples with numeric intensities

Example of raw_data:

| ProteinID | Control\_1 | Control\_2 | Control\_3 | Treatment\_1 | Treatment\_2 | Treatment\_3 |
| --------- | ---------- | ---------- | ---------- | ------------ | ------------ | ------------ |
| Protein1  | 70         | 50         | 27         | 74           | 42           | 57           |
| Protein2  | 0          | 0          | 0          | 55           | 70           | 88           |
| Protein3  | 67         | 52         | 40         | 0            | 0            | 0            |
| Protein4  | 0          | 0          | 0          | 0            | 0            | 0            |


  2. metadata – describes your samples and experimental groups
  Must contains two columns
    - Sample: matches exactly the column names of raw_data (except ProteinID)
    - Group: experimental group/condition

Example of metadata

| Sample       | Group     |
| ------------ | --------- |
| Control\_1   | Control   |
| Control\_2   | Control   |
| Control\_3   | Control   |
| Treatment\_1 | Treatment |
| Treatment\_2 | Treatment |
| Treatment\_3 | Treatment |

# Main Workflow

The recommended workflow follows a logical sequence:
  - proteo.import() – Import and validate raw proteomics data (raw_data) and metadata.
  - proteo.normalize() – Normalize, remove outliers, impute missing values, and filter proteins.
  - proteo.qc() – Generate QC plots (boxplots, PCA) to evaluate consistency.

# Exploratory and differential analysis
  - proteo.corr() – Compute correlations among samples/groups.
  - proteo.dimred() – Perform PCA and UMAP for dimensionality reduction.
  - proteo.enrich() – GO enrichment analysis for protein groups.
  - proteo.heatmap() – Heatmap of top variable proteins.
  - proteo.ttest() – T-test or Mann-Whitney test for protein-level differences.
  - proteo.venn() – Visualize overlaps with Venn diagrams.
  - proteo.volcano() – Volcano plots for differential expression.

## Example (Minimal)

# Load example raw data and metadata
```r
raw_data <- data.frame(
  ProteinID = c("P001","P002","P003"),
  Control1 = c(100,200,150),
  Control2 = c(110,210,145),
  Treatment1 = c(120,180,160),
  Treatment2 = c(115,190,155)
)

metadata <- data.frame(
  Sample = c("Control1","Control2","Treatment1","Treatment2"),
  Group = c("Control","Control","Treatment","Treatment")
)

# Import raw data
obj <- proteo.import(raw_data, metadata) # returns a proteoR object

# Normalize data
normalized_data <- proteo.normalize(raw_data, metadata)

# QC
proteo.qc(normalized_data, metadata)

# Compute correlation
corr_mat <- proteo.corr(normalized_data, metadata)
```

## Contact

For questions, suggestions, or contributions, open an issue or pull request on GitHub.

Thank you for using autoproteoR!
