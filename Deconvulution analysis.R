######DECONVULUTION SVR MANUAL METHOD

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("granulator", update = FALSE, ask = FALSE)

library(granulator)

##  Load Preprocessed Proteomics Data
# Ensure this file has proteins/genes as rows and samples as columns
proteomics_data <- read.csv("proteomics.csv", row.names = 1, check.names = FALSE)

# Also load reference matrix into the environment beforehand! (.csv or seuratObj)
# Example: reference_matrix <- read.csv("reference_matrix.csv", row.names = 1)

## Normalize Expression Data to [0,1]
#-Why Normalize to [0,1]?
#Proteomics data spans a wide dynamic range and varies in scale across both proteins and cell types, normalization is critical.
#Reference matrix (column-wise scaling):
#Each cell type profile is scaled independently. This ensures that no cell type dominates simply because its proteins are measured at higher absolute abundances. It makes the reference signatures comparable across cell types.

#Bulk data (row-wise scaling):
#Each protein is normalized across all bulk samples. This equalizes the contribution of each protein to the regression model, preventing highly abundant proteins (e.g., plasma proteins like albumin or complement factors) from overshadowing lower-abundance, yet cell type–specific proteins.

#Together, this scaling ensures that the regression is driven by relative patterns rather than raw magnitude.

# Column-wise [0,1] scaling for reference matrix (per cell type)

scale_0_1_colwise <- function(x) {
  apply(x, 2, function(col) {
    rng <- range(col, na.rm = TRUE)
    if (diff(rng) == 0 || any(is.na(rng))) return(rep(0, length(col)))
    (col - rng[1]) / diff(rng)
  })
}

# Row-wise scaling for bulk (per gene, across samples)
scale_0_1_rowwise <- function(x) {
  t(apply(x, 1, function(row) {
    rng <- range(row, na.rm = TRUE)
    if (diff(rng) == 0 || any(is.na(rng))) return(rep(0, length(row)))
    (row - rng[1]) / diff(rng)
  }))
}

scaled_proteomics <- scale_0_1_rowwise(proteomics_data)
scaled_reference <- scale_0_1_colwise(reference_matrix)

#Subset to Common Genes
#After normalization, we subset to common genes (or proteins) because deconvolution requires bulk and reference matrices to be defined over the same feature set.
#Without restricting to overlapping features, the regression problem is ill-posed, as the two datasets would be in different spaces. 
#This alignment ensures that only proteins reliably quantified in both bulk and reference contribute to the estimation of cell type proportions.

common_genes <- intersect(rownames(scaled_proteomics), rownames(scaled_reference))
bulk <- scaled_proteomics[common_genes, ]
ref  <- scaled_reference[common_genes, ]

##]: Impute Missing Values (if any)
impute_min <- function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE)
  return(x)
}

bulk_imputed <- t(apply(bulk, 1, impute_min))
ref_imputed  <- t(apply(ref, 1, impute_min))

# Ensure no NAs remain
stopifnot(!anyNA(bulk_imputed), !anyNA(ref_imputed))

# Check dimensions
cat("Genes available:", nrow(ref_imputed), "\nCell types:", ncol(ref_imputed), "\n")

if (nrow(ref_imputed) < ncol(ref_imputed)) {
  stop("ERROR: Number of genes must be >= number of cell types in the reference matrix.")
}


# Check for NA / NaN / Inf in inputs
sum(is.na(bulk_imputed))     # Should be 0
sum(is.na(ref_imputed))      # Should be 0

sum(is.nan(bulk_imputed))    # Should be 0
sum(is.nan(ref_imputed))     # Should be 0

sum(is.infinite(bulk_imputed))  # Should be 0
sum(is.infinite(ref_imputed))   # Should be 0

zero_rows <- rowSums(ref_imputed) == 0
cat("All-zero rows in ref:", sum(zero_rows), "\n")
if (sum(zero_rows) > 0) {
  ref_imputed <- ref_imputed[!zero_rows, ]
  bulk_imputed <- bulk_imputed[!zero_rows, ]
}

apply(ref_imputed, 2, function(col) length(unique(col)))
boxplot(ref_imputed, main = "Distributions of Cell Type Expression in Signature Matrix")

## Run Deconvolution Using SVR ---
result <- deconvolute(
  m = bulk_imputed,
  sigMatrix = ref_imputed,
  methods = "svr",
  use_cores = 1
)

##Inspect
print(head(result$proportions))  # Rows = cell types, Columns = samples

# Save to CSV
write.csv(result$proportions, "cell_type_proportions_granulator.csv")


library(e1071)
#svr manual method 
manual_svr_deconvolve <- function(bulk, ref, nu = 0.5, normalize = TRUE) {
  
  # Check input dimensions
  stopifnot(nrow(bulk) == nrow(ref))
  
  # Initialize result matrix
  proportions <- matrix(NA, nrow = ncol(ref), ncol = ncol(bulk))
  rownames(proportions) <- colnames(ref)
  colnames(proportions) <- colnames(bulk)
  
  # Loop over each sample in bulk
  for (i in seq_len(ncol(bulk))) {
    y <- bulk[, i]
    x <- ref
    
    # Try SVM
    tryCatch({
      model <- svm(x = x, y = y, type = "nu-regression", nu = nu)
      # Approximate cell-type proportions via pseudo-inverse projection
      weights <- coef(lm(y ~ x - 1))  # no intercept
      proportions[, i] <- weights
    }, error = function(e) {
      warning(sprintf("SVM failed for sample %s: %s", colnames(bulk)[i], e$message))
      proportions[, i] <- NA
    })
  }
  
  # Normalize columns to sum to 1
  if (normalize) {
    proportions <- apply(proportions, 2, function(x) {
      if (all(is.na(x))) return(x)
      x / sum(x, na.rm = TRUE)
    })
  }
  
  return(proportions)
}
svr_result <- manual_svr_deconvolve(bulk_filtered, ref_filtered)
head(svr_result)


colnames(bulk_filtered)

bulk_filtered <- bulk_filtered[, !grepl("^X(\\.|$)", colnames(bulk_filtered))]

svr_result <- manual_svr_deconvolve(bulk_filtered, ref_filtered)
head(svr_result)


#plots
# We use pheatmap for heatmaps, ggplot2 for custom plots,
# reshape2 for reshaping matrices, dplyr for summaries,

needed <- c("pheatmap", "ggplot2", "reshape2", "dplyr", "viridis")
for (p in needed) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(pheatmap); library(ggplot2); library(reshape2); library(dplyr); library(viridis)

# svr_result must exist and be a numeric matrix/data.frame
# Rows = cell types, columns = samples
stopifnot("`svr_result` must exist" = exists("svr_result"))
stopifnot("`svr_result` must be a matrix/data.frame" =
            is.matrix(svr_result) || is.data.frame(svr_result))
svr_result <- as.matrix(svr_result)

# Utility functions
# clip_01: force values into [0,1] for display
# normalize_cols: scale each sample to proportions (sum = 1)
# to_long: reshape a matrix into long format for ggplot
clip_01 <- function(x) { x[x < 0] <- 0; x[x > 1] <- 1; x }
normalize_cols <- function(mat) {
  apply(mat, 2, function(col) {
    s <- sum(col, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(col)))
    col / s
  })
}
to_long <- function(mat, value = "Proportion") {
  out <- reshape2::melt(mat, varnames = c("CellType", "Sample"), value.name = value)
  names(out) <- c("CellType","Sample",value)
  out
}

# Define a clean plot theme
# A custom theme for consistency across figures
theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.title = element_text(face = "bold")
    )
}

# Pre-processing
# Clip values for the heatmap, and generate normalized matrices
# - clipped_result: values limited to [0,1]
# - svr_proportion: non-negative, normalized to sum=1 per sample
# - svr_percentage: proportions expressed as %
clipped_result <- clip_01(svr_result)
svr_nonneg     <- pmax(svr_result, 0)
svr_proportion <- normalize_cols(svr_nonneg)
svr_percentage <- svr_proportion * 100

# Assign groups to samples
#samples starting with DR → "DR", PD → "PD", others → "Other"
sample_names <- colnames(svr_result)
sample_group <- ifelse(grepl("^DR", sample_names), "DR",
                       ifelse(grepl("^PD", sample_names), "PD", "Other"))
sample_meta <- data.frame(Sample = sample_names, Group = sample_group)

# Heatmap visualization
# Heatmap of clipped values to highlight global patterns
pheatmap(
  clipped_result,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "steelblue"))(100),
  main = "Cell Type Proportions (Manual SVR) — Clipped [0,1]"
)

# Stacked bar plots
# Convert to long format, merge with sample metadata, then plot
df_prop <- to_long(svr_proportion, value = "Proportion") %>%
  left_join(sample_meta, by = "Sample")
df_pct <- df_prop %>% mutate(Contribution = Proportion * 100)

# (a) Proportions (0–1)
ggplot(df_prop, aes(Sample, Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(
    title = "Cell Type Proportions by Sample",
    subtitle = "Each column normalized to sum to 1",
    x = "Sample", y = "Proportion", fill = "Cell Type"
  ) + theme_pub()

# (b) Percentages (0–100%) faceted by group
ggplot(df_pct, aes(Sample, Contribution, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(
    title = "Cell Type Contributions by Group",
    subtitle = "Faceted by DR / PD / Other",
    x = "Sample", y = "Contribution (%)", fill = "Cell Type"
  ) + theme_pub()

# Group summaries
# Summarize mean ± SD contribution per group & cell type
summary_tbl <- df_pct %>%
  group_by(Group, CellType) %>%
  summarise(
    MeanPct = mean(Contribution, na.rm = TRUE),
    SDPct   = sd(Contribution, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Label = sprintf("%.1f%% ± %.1f", MeanPct, SDPct))
print(summary_tbl)

# Optional: Save outputs
# ggsave("stacked_proportions.png", width = 10, height = 5, dpi = 300)
# ggsave("stacked_percentages_by_group.png", width = 11, height = 5, dpi = 300)
# write.csv(summary_tbl, "group_celltype_summary.csv", row.names = FALSE)

# Session info for reproducibility
# sessionInfo()