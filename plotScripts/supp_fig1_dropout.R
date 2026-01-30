# Calculate and compare gene dropout metrics across datasets.

library(dplyr)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Load Data
# --------------------------------------------------------------------------------
# Load probe genes
genes_143 <- read.csv(file.path(input_dir, "Mouse_panel_143_gene_R_L_final.csv"))
genes_143$gene_short_name <- sapply(strsplit(genes_143$probe_id, "_"), function(x) x[2])
# head(genes_143) # Optional check
probe_genes <- unique(genes_143$gene_id) 

# Load CDS objects
cds_on <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-2-TRTL-ON_processed_cds_object.rds"))
cds_dt <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_16krpc_cds_object.rds"))
cds_dt_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_7krpc_cds_object.rds"))
cds_TRTL_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_7krpc_cds_object.rds"))

cds_list <- list(
  "cds_on" = cds_on,
  "cds_dt" = cds_dt,
  "cds_dt_7000" = cds_dt_7000rpc,
  "cds_TRTL_7000" = cds_TRTL_7000rpc
)

# --------------------------------------------------------------------------------
# Calculate Dropout Metrics
# --------------------------------------------------------------------------------
results <- lapply(names(cds_list), function(name) {
  cds <- cds_list[[name]]
  
  # 1. Access the RAW COUNT matrix (integers)
  raw_counts <- counts(cds)
  
  # 2. Subset to only the probe genes
  # We match the probe_genes list to the 'id' column in rowData
  if(!exists("probe_genes")) stop("probe_genes variable is not defined.")
  
  probe_indices <- which(rowData(cds)$id %in% probe_genes)
  sub_counts <- raw_counts[probe_indices, ]
  
  # 3. Calculate metrics based on RAW integers
  # Cells per gene with > 1 count (High-confidence detection)
  cells_with_gt1 <- Matrix::rowSums(sub_counts > 1)
  
  # Cells per gene with >= 1 count (Any detection)
  cells_with_ge1 <- Matrix::rowSums(sub_counts >= 1)
  
  num_cells_total <- ncol(cds)
  
  # 4. Aggregate results
  # We calculate how many genes meet your 0.5% frequency threshold 
  # based on the "Strict" (>1 count) requirement.
  df_metrics <- data.frame(
    gene_id = rowData(cds)$id[probe_indices],
    strict_detect = cells_with_gt1
  ) %>%
    mutate(pct_strict = (strict_detect / num_cells_total) * 100)
  
  data.frame(
    dataset = name,
    `Total Cells` = num_cells_total,
    `Genes Detected (> 1 Count)` = sum(df_metrics$strict_detect >= 1),
    `Genes Detected (> 1 Count) in at least 0.5% of total cells` = sum(df_metrics$pct_strict > 0.5),
    check.names = FALSE
  )
})

comparison_table <- do.call(rbind, results)
print(comparison_table)

write.csv(comparison_table, 
          file.path(output_dir, "dropout_comparison_table.csv"), 
          row.names = FALSE)

