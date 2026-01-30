suppressPackageStartupMessages({
  library(tidyr)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(monocle3)
  library(gridExtra)
  library(ggrepel)
  library(data.table)  
  library(pheatmap)
  library(patchwork)
  library(SummarizedExperiment)
})

# Reading in arguments from terminal
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide input_dir and output_dir as arguments.")
}
input_dir <- args[1]
output_dir <- args[2]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

allen_dir <- file.path(output_dir, "allen")
if (!dir.exists(allen_dir)) {
  dir.create(allen_dir, recursive = TRUE)
}

# --- Load Data ---

# 16k DT (Full Dataset)
# cds_dt from fig1_major_cell_types 
 
cds_dt_path <- file.path(input_dir, "cds/cds_dt_cell_type_annotated/cds_object.rds")
if (!file.exists(cds_dt_path)) {
  # Fallback to aligned UMAP if annotated not found
  cds_dt_path <- file.path(input_dir, "cds/cds_dt_aligned_UMAP/cds_object.rds") 
}
if (!file.exists(cds_dt_path)) stop("Could not find cds_dt input.")

cds_dt <- readRDS(cds_dt_path)

# 7k RPC DT (Subsampled)
# From fig1_7k_downsampled.R -> output/cds/cds_dt_7krpc_aligned_UMAP
cds_dt_7krpc_path <- file.path(input_dir, "cds/cds_dt_7krpc_aligned_UMAP/cds_object.rds")
if (!file.exists(cds_dt_7krpc_path)) stop("Could not find cds_dt_7krpc input.")
cds_dt_7krpc <- readRDS(cds_dt_7krpc_path)

# 7k RPC FF (Subsampled)
# From fig1_7k_downsampled.R -> output/cds/cds_ff_7krpc_aligned_UMAP
cds_ff_7krpc_path <- file.path(input_dir, "cds/cds_ff_7krpc_aligned_UMAP/cds_object.rds")
if (!file.exists(cds_ff_7krpc_path)) stop("Could not find cds_ff_7krpc input.")
cds_ff_7krpc <- readRDS(cds_ff_7krpc_path)


# --- Define Save Function ---

save_count_matrix <- function(cds, output_path, filename = "count_matrix.csv") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("install the 'data.table' package.")
  }
  
  # Extract and transpose count matrix
  t_counts <- as.matrix(exprs(cds))
  count_matrix <- t(t_counts)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # Write to CSV
  out_file <- file.path(output_path, filename)
  fwrite(as.data.frame(count_matrix), file = out_file, quote = FALSE, row.names = TRUE)
  
  message("Count matrix saved to: ", out_file)
}

# --- Export Matrices ---

# FF 7k
message("Saving FF 7k matrix...")
save_count_matrix(cds = cds_ff_7krpc, output_path = allen_dir, filename = "count_matrix_ff.csv")

# Code to split large FF matrix (if needed, consistent with reference script logic)
ff_file <- file.path(allen_dir, "count_matrix_ff.csv")
if (file.exists(ff_file) && file.size(ff_file) > 2e9) { # Check if > 2GB (approx) 
  message("Splitting FF matrix...")
  full_data <- fread(ff_file, data.table = FALSE)
  rownames(full_data) <- full_data[[1]]
  full_data <- full_data[, -1]
  
  num_rows <- nrow(full_data)
  mid_point <- floor(num_rows / 2)
  
  part1 <- full_data[1:mid_point, ]
  part2 <- full_data[(mid_point + 1):num_rows, ]
  
  fwrite(part1, file.path(allen_dir, "count_matrix_ff_part1.csv"), quote = FALSE, row.names = TRUE)
  fwrite(part2, file.path(allen_dir, "count_matrix_ff_part2.csv"), quote = FALSE, row.names = TRUE)
}

# DT 7k
message("Saving DT 7k matrix...")
save_count_matrix(cds = cds_dt_7krpc, output_path = allen_dir, filename = "count_matrix_dt_7k.csv")

# DT 16k (Main)
message("Saving DT 16k matrix...")
save_count_matrix(cds = cds_dt, output_path = allen_dir, filename = "count_matrix_dt_16k.csv")


# --- Merge Logic (Placeholder/Template) ---

# This section assumes the user has run the Allen mapping externally and placed the results back in the output folder.

merge_allen_metadata <- function(cds, annotation_csv_path, skip_lines = 4) {
  if (!file.exists(annotation_csv_path)) {
    warning(paste("Annotation file not found:", annotation_csv_path))
    return(cds)
  }
  
  # Read in Allen Brain Map annotations
  map_allen <- read_csv(annotation_csv_path, skip = skip_lines)
  
  # Rename column for merging
  if ("cell_id" %in% names(map_allen)) {
      map_allen <- map_allen %>% dplyr::rename(cell = cell_id)
  }
  
  # Prepare colData for merging
  colData(cds)$cell <- rownames(colData(cds))
  
  # Merge annotation with colData
  merged_data <- left_join(as.data.frame(colData(cds)), map_allen, by = "cell")
  rownames(merged_data) <- merged_data$cell
  
  # Replace colData with merged data
  colData(cds) <- DataFrame(merged_data)
  
  return(cds)
}

# --- Merge Mapped Results ---

# FF 7k: Combine split mapped files if they exist, or read combined
mapped_ff_path <- file.path(allen_dir, "mapped_ff_combined.csv")

# Logic to combine parts if main doesn't exist but parts do (Allen output usually has specific names)
# We look for files containing "count_matrix_ff_part" and ".csv"
part1_files <- list.files(allen_dir, pattern = "count_matrix_ff_part1.*\\.csv", full.names = TRUE)
part2_files <- list.files(allen_dir, pattern = "count_matrix_ff_part2.*\\.csv", full.names = TRUE)

if (!file.exists(mapped_ff_path) && length(part1_files) > 0 && length(part2_files) > 0) {
    message("Combining mapped FF parts...")
    # Assuming the first match is correct if multiple exist
    mapped1 <- fread(part1_files[1], skip = 4, fill = TRUE, header = TRUE) 
    mapped2 <- fread(part2_files[1], skip = 4, fill = TRUE, header = TRUE)
    mapped_combined <- rbind(mapped1, mapped2)
    fwrite(mapped_combined, mapped_ff_path, quote = FALSE, row.names = FALSE)
}

message("Merging Allen metadata for FF 7k...")
cds_ff_7krpc <- merge_allen_metadata(
  cds = cds_ff_7krpc, 
  annotation_csv_path = mapped_ff_path, 
  skip_lines = 0
)

# DT 7k
mapped_dt_7k_path <- file.path(allen_dir, "mapped_dt_7k.csv")
# Look for original download name if generic not found
mapped_dt_7k_files <- list.files(allen_dir, pattern = "count_matrix_dt_7k.*\\.csv", full.names = TRUE)
if (!file.exists(mapped_dt_7k_path) && length(mapped_dt_7k_files) > 0) {
    mapped_dt_7k_path <- mapped_dt_7k_files[1]
}

message("Merging Allen metadata for DT 7k...")
cds_dt_7krpc <- merge_allen_metadata(
  cds = cds_dt_7krpc,
  annotation_csv_path = mapped_dt_7k_path,
  skip_lines = 4 # Default for Allen downloads
)

# DT 16k
mapped_dt_16k_path <- file.path(allen_dir, "mapped_dt_16k.csv")
mapped_dt_16k_files <- list.files(allen_dir, pattern = "count_matrix_dt_16k.*\\.csv", full.names = TRUE)
if (!file.exists(mapped_dt_16k_path) && length(mapped_dt_16k_files) > 0) {
    mapped_dt_16k_path <- mapped_dt_16k_files[1]
}

message("Merging Allen metadata for DT 16k...")
cds_dt <- merge_allen_metadata(
  cds = cds_dt,
  annotation_csv_path = mapped_dt_16k_path,
  skip_lines = 4
)
# class name grouping using helper script
dt_grouped<-group_cell_subtypes(as.data.frame(colData(cds_dt)))
dt_7krpc_grouped<-group_cell_subtypes(as.data.frame(colData(cds_dt_7krpc)))
ff_7krpc_grouped<-group_cell_subtypes(as.data.frame(colData(cds_ff_7krpc)))

colData(cds_dt)$major_group_cell_type<-dt_grouped$major_group_cell_type
colData(cds_dt_7krpc)$major_group_cell_type<-dt_7krpc_grouped$major_group_cell_type
colData(cds_ff_7krpc)$major_group_cell_type<-ff_7krpc_grouped$major_group_cell_type

# Save final objects
if (!dir.exists(file.path(output_dir, "cds"))) dir.create(file.path(output_dir, "cds"), recursive = TRUE)

saveRDS(cds_dt_7krpc, file.path(output_dir, "cds/sci-RNA-TRTL-mouse-1-2-combined_processed_dT_7krpc_cds_object.rds"))
saveRDS(cds_ff_7krpc, file.path(output_dir, "cds/sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_7krpc_cds_object.rds"))
saveRDS(cds_dt, file.path(output_dir, "cds/sci-RNA-TRTL-mouse-1-2-combined_processed_dT_16krpc_cds_object.rds"))

