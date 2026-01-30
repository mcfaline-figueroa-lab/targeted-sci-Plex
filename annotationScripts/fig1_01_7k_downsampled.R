# Note: The individual subsampled cds's per condition are not deposited in GEO as there were so many, but can be provided upon request

suppressPackageStartupMessages({
  library(tidyr)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(monocle3)
  library(gridExtra)
  library(ggrepel)
  library(Seurat)
  library(SeuratObject)
  library(RColorBrewer)
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

if (!dir.exists(file.path(output_dir, "UMAPs"))) {
  dir.create(file.path(output_dir, "UMAPs"), recursive = TRUE)
}

# --- 1. Data Assembly downsampled 7k rpc---

cds.pre.path.1.dt = file.path(input_dir, "cds_precell_prehash_dt_exp1.RDS")
cds.pre.path.1.targ = file.path(input_dir, "cds_precell_prehash_targ_exp1.RDS")
cds.pre.path.2.dt = file.path(input_dir, "cds_precell_prehash_dt_exp2.RDS")
cds.pre.path.2.targ = file.path(input_dir, "cds_precell_prehash_ff_exp2.RDS") # Renamed from '55' for consistency
cds.pre.path.on = file.path(input_dir, "cds_precell_prehash_ON.RDS")
# ff refers to TRTL at 55C 7 min anneal as opposed to TRTL-ON

cds.pre.dt.1 <- readRDS(cds.pre.path.1.dt)
cds.pre.ff.1 <- readRDS(cds.pre.path.1.targ)
cds.pre.dt.2 <- readRDS(cds.pre.path.2.dt)
cds.pre.ff.2 <- readRDS(cds.pre.path.2.targ)
cds.pre.on <- readRDS(cds.pre.path.on)

colData(cds.pre.dt.1)$method <- "dt"; colData(cds.pre.dt.1)$experiment <- "1"
colData(cds.pre.ff.1)$method <- "ff"; colData(cds.pre.ff.1)$experiment <- "1"
colData(cds.pre.dt.2)$method <- "dt"; colData(cds.pre.dt.2)$experiment <- "2"
colData(cds.pre.ff.2)$method <- "ff"; colData(cds.pre.ff.2)$experiment <- "2"
colData(cds.pre.on)$method <- "on"; colData(cds.pre.on)$experiment <- "2"

cds_ff <- combine_cds(list(cds.pre.ff.1, cds.pre.ff.2), cell_names_unique = TRUE)
cds_dt <- combine_cds(list(cds.pre.dt.1, cds.pre.dt.2), cell_names_unique = TRUE)

names(colData(cds_ff))[names(colData(cds_ff)) == "Cell"] <- "cell"
names(colData(cds_dt))[names(colData(cds_dt)) == "Cell"] <- "cell"
names(colData(cds.pre.on))[names(colData(cds.pre.on)) == "Cell"] <- "cell"

# Replicate Metadata
RT_samplesheet <- read_csv(file.path(input_dir, "replicate_mouse_exp2_250825_all_3.csv"))

process_metadata <- function(cds, samplesheet) {
  cols <- colData(cds) %>% 
    as.data.frame() %>%
    mutate(cell_id = cell) %>%
    tidyr::separate(col = cell_id, into = c('P7', 'P5', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
    mutate(RT = paste('RT_BC_', n1, sep = '')) %>%
    mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
    dplyr::select(-n1, -n2)
  
  cols <- cols %>%
    left_join(samplesheet, by = c("RT" = "RT_cell_map")) %>%
    mutate(replicate = ifelse(experiment == "1", "1", replicate)) # Adjust replicate for experiment 1
  
  return(cols)
}

cols_ff <- process_metadata(cds_ff, RT_samplesheet)
cols_dt <- process_metadata(cds_dt, RT_samplesheet)
# ON data logic
cols_on <- colData(cds.pre.on) %>% 
  as.data.frame() %>%
  mutate(cell_id = cell) %>%
  tidyr::separate(col = cell_id, into = c('P7', 'P5', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
  mutate(RT = paste('RT_BC_', n1, sep = '')) %>%
  mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
  select(-n1, -n2) %>%
  left_join(RT_samplesheet, by = c("RT" = "RT_cell_map"))

# Process CDS with Columns
process_cds_with_cols <- function(cols, cds) {
  rownames(cols) <- cols$cell
  cols <- DataFrame(cols)
  
  if (!identical(rownames(cols), rownames(colData(cds)))) {
    stop("Rownames of `cols` and `cds` colData do not match.")
  }
  
  colData(cds) <- cols
  colData(cds)$log10.umi <- log10(colData(cds)$n.umi)
  
  mt_genes <- rowData(cds) %>% as.data.frame() %>% filter(grepl("^MT-", gene_short_name))
  mt_gene_ids <- rownames(mt_genes)
  colData(cds)$percent_mito <- 100 * (colSums(exprs(cds)[mt_gene_ids, , drop = FALSE]) / colSums(exprs(cds)))
  
  return(cds)
}

cds_ff <- process_cds_with_cols(cols_ff, cds_ff)
cds_dt <- process_cds_with_cols(cols_dt, cds_dt)
cds_on <- process_cds_with_cols(cols_on, cds.pre.on)

# --- QC & Filtering ---

analyze_cds_quality <- function(cds, plot_knee = TRUE, umi_cutoff = NULL, plot_title = NULL, output_path = NULL) {
  cds <- detect_genes(cds)
  cds <- estimate_size_factors(cds)
  cds <- detect_genes(cds)

  if (plot_knee && !is.null(output_path)) {
    knee_df <- as.data.frame(colData(cds)) %>%
      mutate(cells = seq_along(cell), umis = n.umi) %>%
      filter(!is.na(replicate)) %>%    
      select(cells, umis, replicate)
    
    knee_df_sorted <- knee_df %>%
      group_by(replicate) %>%
      arrange(replicate, desc(umis)) %>%
      mutate(cells_sorted = row_number()) %>%
      ungroup()
    
    knee_plot <- ggplot(knee_df_sorted, aes(x = log10(cells_sorted), y = log10(umis))) +
      geom_point(alpha = 0.5) +
      facet_wrap(~ replicate) +
      theme_minimal() +
      labs(x = "Log10(Cells)", y = "Log10(UMIs)") +
      theme(legend.position = "none")
    
    if (!is.null(umi_cutoff)) {
      knee_plot <- knee_plot + 
        geom_hline(yintercept = log10(umi_cutoff), linetype = "dashed", color = "red") +
        annotate("text", x = 0, y = log10(umi_cutoff), label = paste0(umi_cutoff, " UMIs"), 
                 hjust = -0.1, vjust = -0.5, color = "red", size = 3)
    }
    if (!is.null(plot_title)) {
      knee_plot <- knee_plot + ggtitle(plot_title)
    }
    
    if (!dir.exists(dirname(output_path))) dir.create(dirname(output_path), recursive = TRUE)
    ggsave(output_path, plot = knee_plot, height = 2, width = 7, dpi = 600)
  }
  return(cds)
}

cds_dt <- analyze_cds_quality(cds_dt, umi_cutoff=700, plot_title = "dT randomN RT", output_path = file.path(output_dir, "plots/dt_knees_subsampled.png"))
cds_ff <- analyze_cds_quality(cds_ff, umi_cutoff=300, plot_title ="targeted RT", output_path = file.path(output_dir, "plots/ff_knees_subsampled.png"))
cds_on <- analyze_cds_quality(cds_on, umi_cutoff=70, plot_title ="targeted ON RT", output_path = file.path(output_dir, "plots/on_knees_subsampled.png"))

# Combined Knee Plot Code
umi_cutoffs <- data.frame(
  cds_source = c("cds_on", "cds_dt", "cds_ff"),
  cutoff = c(70, 700, 300)
)

combined_knee_df <- bind_rows(
  as.data.frame(colData(cds_on)) %>% mutate(umis = n.umi, cds_source = "cds_on") %>% select(umis, cds_source),
  as.data.frame(colData(cds_dt)) %>% mutate(umis = n.umi, cds_source = "cds_dt") %>% select(umis, cds_source),
  as.data.frame(colData(cds_ff)) %>% mutate(umis = n.umi, cds_source = "cds_ff") %>% select(umis, cds_source)
) %>%
  group_by(cds_source) %>%
  arrange(desc(umis)) %>%
  mutate(cells_sorted = row_number()) %>%
  ungroup() %>%
  left_join(umi_cutoffs, by = "cds_source")

annotation_df <- combined_knee_df %>% distinct(cds_source, cutoff)
cds_labels <- c("cds_dt" = "DT", "cds_ff" = "FF", "cds_on" = "ON")

combined_knee_plot <- ggplot(combined_knee_df, aes(x = log10(cells_sorted), y = log10(umis))) +
  geom_point(alpha = 0.5) +
  geom_hline(aes(yintercept = log10(cutoff)), linetype = "dashed", color = "red") +
  facet_wrap(~ cds_source, labeller = as_labeller(cds_labels)) +
  theme_minimal() +
  geom_text(data = annotation_df, aes(x = 0, y = log10(cutoff), label = paste0(cutoff, " UMIs")), hjust = -0.1, vjust = -0.5, color = "red", size = 3)

ggsave(file.path(output_dir, "plots/all_method_knees.png"), plot = combined_knee_plot, height=2, width=5, dpi = 900)

# --- Preprocessing with Probe Genes ---

genes_143 <- read.csv(file.path(input_dir, "Mouse_panel_143_gene_R_L_final.csv"))
genes_143$gene_short_name <- sapply(strsplit(genes_143$probe_id, "_"), function(x) x[2])
probe_genes <- unique(genes_143$gene_id)

filter_and_preprocess_cds <- function(cds, umi_cutoff = NULL, probe_genes = NULL, save_plot_path = NULL) {
  cds <- cds[, colData(cds)$n.umi > umi_cutoff]
  
  expressed_genes <- rowData(cds) %>% 
    as.data.frame() %>% 
    filter(num_cells_expressed >= 0.05 * ncol(cds)) %>% 
    pull(id) %>% 
    unique()

  if (!is.null(probe_genes)) {
    use_genes <- unique(c(probe_genes, expressed_genes))
  } else {
    use_genes <- expressed_genes
  }
  
  cds <- preprocess_cds(cds, num_dim = 30, use_genes = use_genes)
  
  if (!is.null(save_plot_path)) {
    p <- plot_pc_variance_explained(cds) + geom_vline(xintercept = 15) + theme_minimal()
    if (!dir.exists(dirname(save_plot_path))) dir.create(dirname(save_plot_path), recursive = TRUE)
    ggsave(save_plot_path, plot = p, width = 2, height = 2, dpi = 900)
  }
  return(cds)
}

cds_ff_7krpc <- filter_and_preprocess_cds(cds_ff, umi_cutoff = 300, probe_genes = probe_genes, save_plot_path = file.path(output_dir, "QC/elbow_ff.png"))
cds_dt_7krpc <- filter_and_preprocess_cds(cds_dt, umi_cutoff = 700, probe_genes = probe_genes, save_plot_path = file.path(output_dir, "QC/elbow_dt.png"))

# Reductions & Alignment
# DT
cds_dt_7krpc <- detect_genes(cds_dt_7krpc)
cds_dt_7krpc <- preprocess_cds(cds_dt_7krpc, num_dim = 15, use_genes = c(probe_genes)) # Just probe genes per original script section (or combined if preferred)
cds_dt_7krpc <- reduce_dimension(cds_dt_7krpc, preprocess_method = "PCA", umap.min_dist = 0.1, umap.n_neighbors = 25)
cds_dt_7krpc <- align_cds(cds_dt_7krpc, preprocess_method = "PCA", alignment_group ="experiment")
cds_dt_7krpc <- reduce_dimension(cds_dt_7krpc, preprocess_method = "Aligned", umap.min_dist = 0.1, umap.n_neighbors = 25)
p_dt <- plot_cells(cds_dt_7krpc, scale_to_range = FALSE, graph_label_size = 20, color_cells_by = "experiment") + facet_wrap(~experiment)  
ggsave(file.path(output_dir, "QC/aligned_UMAP_dt.png"), plot = p_dt, height = 5, width = 6, dpi=600)

# FF
cds_ff_7krpc <- preprocess_cds(cds_ff_7krpc, num_dim = 15, use_genes = c(probe_genes))
cds_ff_7krpc <- reduce_dimension(cds_ff_7krpc, preprocess_method = "PCA", umap.min_dist = 0.1, umap.n_neighbors = 25)
cds_ff_7krpc <- align_cds(cds_ff_7krpc, preprocess_method = "PCA", alignment_group ="experiment")
cds_ff_7krpc <- reduce_dimension(cds_ff_7krpc, preprocess_method = "Aligned", umap.min_dist = 0.1, umap.n_neighbors = 25)
p_ff <- plot_cells(cds_ff_7krpc, scale_to_range = FALSE, graph_label_size = 20, color_cells_by = "experiment") + facet_wrap(~experiment)
ggsave(file.path(output_dir, "QC/aligned_UMAP_targeted.png"), plot = p_ff, height = 5, width = 6, dpi=600)

save_monocle_objects(cds_dt_7krpc, file.path(output_dir, "cds/cds_dt_7krpc_aligned_UMAP"))
save_monocle_objects(cds_ff_7krpc, file.path(output_dir, "cds/cds_ff_7krpc_aligned_UMAP"))



