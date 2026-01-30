suppressPackageStartupMessages({
  library(glmGamPoi)
  library(tidyr)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(gridExtra)
  library(ggrepel)
  library(Seurat)
  library(SeuratObject)
  library(monocle3)
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

# --- 1. Data Assembly ---

cds.pre.path.1.dt = file.path(input_dir, "sci-RNA-TRTL-mouse-1_cds_precell_prehash.RDS")
cds.pre.path.2.dt = file.path(input_dir, "sci-RNA-TRTL-mouse-2_cds_precell_prehash.RDS")

cds.dt.1 <- readRDS(cds.pre.path.1.dt)
cds.dt.2 <- readRDS(cds.pre.path.2.dt)
#not yet only dt, have to filter after getting p7 into metadata below:

names(colData(cds.dt.1))[names(colData(cds.dt.1)) == "Cell"] <- "cell"
names(colData(cds.dt.2))[names(colData(cds.dt.2)) == "Cell"] <- "cell"
colData(cds.dt.1)$experiment <- "1"
colData(cds.dt.2)$experiment <- "2"
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
    mutate(replicate = ifelse(experiment == "1", "1", replicate))
  
  return(cols)
}
cols.dt.1 <- process_metadata(cds.dt.1, RT_samplesheet)
cols.dt.2 <- process_metadata(cds.dt.2, RT_samplesheet)

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

cds.dt.1 <- process_cds_with_cols(cols.dt.1, cds.dt.1)
cds.dt.2 <- process_cds_with_cols(cols.dt.2, cds.dt.2)

#assign dt and TRTL from P7 and filter for just dt
colData(cds.dt.1) <- DataFrame(as.data.frame(colData(cds.dt.1)) %>%
  mutate(method = ifelse(P7 %in% c("01A", "02A", "03A", "04A", "05A", "06A"), "TRTL", "dt")))
cds.dt.1 <- cds.dt.1[, colData(cds.dt.1)$method == "dt"]

colData(cds.dt.2) <- DataFrame(as.data.frame(colData(cds.dt.2)) %>%
  mutate(method = ifelse(P7 %in% c("01C", "02C", "03C", "04C", "05C", "06C", "07C", "08C"), "TRTL", "dt")))
cds.dt.2 <- cds.dt.2[, colData(cds.dt.2)$method == "dt"]
cds_dt <- combine_cds(list(cds.dt.1, cds.dt.2), cell_names_unique = TRUE)

colData(cds_dt)$log10.umi <- log10(colData(cds_dt)$n.umi)

# --- 2. Preprocessing & Alignment ---

expressed_genes <- rowData(cds_dt) %>% 
  as.data.frame() %>% 
  filter(num_cells_expressed >= 0.05 * ncol(cds_dt)) %>% 
  pull(id) %>% 
  unique()

cds_dt <- preprocess_cds(cds_dt,
                         num_dim = 15,
                         use_genes = c(expressed_genes)) 

cds_dt <- align_cds(cds_dt,
                    preprocess_method = "PCA",
                    alignment_group ="experiment")

cds_dt <- reduce_dimension(cds_dt,
                           preprocess_method = "Aligned",
                           umap.min_dist = 0.1,
                           umap.n_neighbors = 20,
                           umap.metric = "cosine")

plot_cells(cds_dt, color_cells_by = "experiment")
ggsave(file.path(output_dir, "UMAPs/experiments_aligned.png"), height = 4, width = 4, dpi = 900)

# --- 3. Clustering ---

cds_dt <- cluster_cells(cds_dt, reduction_method = "Aligned", cluster_method = "leiden", resolution = 2e-3) 
colData(cds_dt)$Cluster <- monocle3::clusters(cds_dt, reduction_method="Aligned")

plot_cells(cds_dt, color_cells_by = "Cluster", cell_size = 0.5, group_label_size = 3) + theme_void() 
ggsave(file.path(output_dir, "UMAPs/leiden_UMAP_clusters_Aligned.png"), height = 4, width = 4, dpi = 900)

# --- 4. Marker Analysis ---

marker_lists <- list(
  Astrocyte = c("Gfap","Aqp4","S100b"),
  Neuron = c("Rbfox3","Mapt","Snap25","Car10","Kcnd2"),
  Microglia = c("Trem2","Tmem119"),
  Endothelial = c("Cldn5","Emcn","Esam","Flt1"),
  Oligodendrocyte = c("Cnp","Olig1","Olig2","Sox10")
)

gene_symbols <- rowData(cds_dt)$gene_short_name
gene_ids <- rownames(cds_dt)
symbol_to_id <- setNames(gene_ids, gene_symbols)

marker_lists_ids <- lapply(marker_lists, function(x) symbol_to_id[x])
marker_lists_ids <- lapply(marker_lists_ids, function(x) x[!is.na(x)])
all_ids <- unique(unlist(marker_lists_ids))

norm_mat <- normalized_counts(cds_dt)
norm_markers <- norm_mat[all_ids, ]

scaled <- t(scale(t(norm_markers)))
scaled[scaled > 2] <- 2
scaled[scaled < -2] <- -2
scaled[is.na(scaled)] <- 0

gene_order <- unlist(marker_lists_ids)
gene_order <- gene_order[gene_order %in% rownames(scaled)]
scaled <- scaled[gene_order, ]

cd <- colData(cds_dt)
set.seed(10)

# Marker heatmap
cells_by_cluster <- split(colnames(cds_dt), cd$Cluster)
subsampled_cells <- unlist(
  lapply(cells_by_cluster, function(cells) {
    if (length(cells) > 200) sample(cells, 200) else cells
  }),
  use.names = FALSE
)

cluster_order <- cd$Cluster[subsampled_cells]
cell_order <- order(cluster_order)
final_cells <- subsampled_cells[cell_order]
scaled <- scaled[, final_cells]

annotation_col <- data.frame(Cluster = cd$Cluster[final_cells])
rownames(annotation_col) <- final_cells

blue_red_palette <- colorRampPalette(c("blue", "white", "red"))(50)
marker_groups <- names(marker_lists)
num_groups <- length(marker_groups)
marker_colors <- setNames(brewer.pal(max(3, min(8, num_groups)), "Dark2")[seq_len(num_groups)], marker_groups)
annotation_colors <- list(marker_group = marker_colors)

pheatmap(
  scaled,
  color = blue_red_palette,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row = 6,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  filename = file.path(output_dir, "UMAPs/marker_heatmap.png"),
  width = 8,
  height = 10,
  dpi = 900
)

# --- 5. Cell Type Assignment ---

colData(cds_dt)$Cluster <- as.character(colData(cds_dt)$Cluster)

colData(cds_dt)$major_cell_types <- case_when(
  colData(cds_dt)$Cluster %in% c("10","12") ~ "Astrocyte",
  colData(cds_dt)$Cluster == "13" ~ "Endothelial",
  colData(cds_dt)$Cluster %in% c("6") ~ "Oligodendrocyte",
  colData(cds_dt)$Cluster == "16" ~ "Microglia",
  TRUE ~ "Neuron"
)

# --- 6. Seurat Integration & Label Transfer ---

# Prepare Reference (so_dt)
library(future)
plan(sequential)

exprs_cds_dt <- exprs(cds_dt)
metadata_cds_dt <- colData(cds_dt)
so_dt <- CreateSeuratObject(exprs_cds_dt, meta.data = metadata_cds_dt %>% as.data.frame(), assay = "RNA")
so_dt <- NormalizeData(so_dt)
so_dt$experiment <- as.factor(so_dt$experiment)

so_dt[["RNA"]] <- split(so_dt[["RNA"]], f = so_dt$experiment)
so_dt <- SCTransform(so_dt)
so_dt <- RunPCA(so_dt, assay = "SCT")
so_dt <- RunUMAP(so_dt, dims = 1:20)
so_dt <- IntegrateLayers(object = so_dt, method = JointPCAIntegration, normalization.method = "SCT", verbose = F)
so_dt <- FindNeighbors(so_dt, reduction = "integrated.dr", dims = 1:20)
so_dt <- RunUMAP(so_dt, dims = 1:20, reduction = "integrated.dr")
saveRDS(so_dt, file.path(output_dir, "so_dt.rds"))

# --- Target 1: DT 7000 RPC ---

cds_dt_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_7krpc_cds_object.rds"))
# Note: this already has the cell type annotation if you just want to run this code and verify annotations, but use the output of fig1_7krpc for the preannotation cds

exprs_cds_dt_7000rpc <- exprs(cds_dt_7000rpc)
metadata_cds_dt_7000rpc <- colData(cds_dt_7000rpc)
so_dt_7000rpc <- CreateSeuratObject(exprs_cds_dt_7000rpc, meta.data = metadata_cds_dt_7000rpc %>% as.data.frame(), assay = "RNA")
so_dt_7000rpc <- NormalizeData(so_dt_7000rpc)
so_dt_7000rpc <- FindVariableFeatures(so_dt_7000rpc, selection.method = "vst", nfeatures = 1000)
so_dt_7000rpc[["RNA"]] <- split(so_dt_7000rpc[["RNA"]], f = so_dt_7000rpc$experiment)
so_dt_7000rpc <- SCTransform(so_dt_7000rpc)
so_dt_7000rpc <- RunPCA(so_dt_7000rpc, assay = "SCT")
so_dt_7000rpc <- RunUMAP(so_dt_7000rpc, dims = 1:20)
so_dt_7000rpc <- IntegrateLayers(object = so_dt_7000rpc, method = JointPCAIntegration, normalization.method = "SCT", verbose = F)
so_dt_7000rpc <- FindNeighbors(so_dt_7000rpc, reduction = "integrated.dr", dims = 1:20)
so_dt_7000rpc <- RunUMAP(so_dt_7000rpc, dims = 1:20, reduction = "integrated.dr")

so.anchors <- FindTransferAnchors(
    reference = so_dt,
    query = so_dt_7000rpc,
    dims = 1:25,
    reference.assay = "RNA", 
    query.assay = "RNA",
    features = unique(c(probe_genes, VariableFeatures(so_dt, assay = "SCT"))), 
    reference.reduction = "integrated.dr" 
)
predictions <- TransferData(anchorset = so.anchors, refdata = so_dt$major_cell_types, dims = 1:15)
so_dt_7000rpc <- AddMetaData(so_dt_7000rpc, metadata = predictions)

p1 <- DimPlot(so_dt_7000rpc, group.by = "predicted.id") + facet_wrap("predicted.id")
ggsave(file.path(output_dir, "UMAPs/so_dt_7000rpc_predicted.png"), plot = p1, height = 4, width = 6, dpi = 300)
saveRDS(so_dt_7000rpc, file.path(output_dir, "so_dt_7000rpc_integrated_experiment_2.rds"))


# --- Target 2: TRTL 7000 RPC (FF) ---

cds_ff_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_7krpc_cds_object.rds"))
# Note: this already has the cell type annotation if you just want to run this code and verify annotations, but use the output of fig1_7krpc for the preannotation cds

exprs_cds_ff_7000rpc <- exprs(cds_ff_7000rpc)
metadata_cds_ff_7000rpc <- colData(cds_ff_7000rpc)
so_ff_7000rpc <- CreateSeuratObject(exprs_cds_ff_7000rpc, meta.data = metadata_cds_ff_7000rpc %>% as.data.frame(), assay = "RNA")

so_ff_7000rpc[["RNA"]] <- split(so_ff_7000rpc[["RNA"]], f = so_ff_7000rpc$experiment)
so_ff_7000rpc <- SCTransform(so_ff_7000rpc) 
so_ff_7000rpc <- RunPCA(so_ff_7000rpc, assay = "SCT")
so_ff_7000rpc <- RunUMAP(so_ff_7000rpc, dims = 1:20)
so_ff_7000rpc <- IntegrateLayers(object = so_ff_7000rpc, method = JointPCAIntegration, normalization.method = "SCT", verbose = F)
so_ff_7000rpc <- FindNeighbors(so_ff_7000rpc, reduction = "integrated.dr", dims = 1:20)
so_ff_7000rpc <- RunUMAP(so_ff_7000rpc, dims = 1:20, reduction = "integrated.dr")

so.anchors <- FindTransferAnchors(
    reference = so_dt,
    query = so_ff_7000rpc,
    dims = 1:25,
    reference.assay = "RNA", 
    query.assay = "RNA",
    features = unique(c(probe_genes, VariableFeatures(so_dt, assay = "SCT")[1:700])),
    reference.reduction = "integrated.dr"
)
predictions <- TransferData(anchorset = so.anchors, refdata = so_dt$major_cell_types, dims = 1:15)
so_ff_7000rpc <- AddMetaData(so_ff_7000rpc, metadata = predictions)

p2 <- DimPlot(so_ff_7000rpc, group.by = "predicted.id") + facet_wrap("predicted.id")
ggsave(file.path(output_dir, "UMAPs/so_ff_7000rpc_predicted.png"), plot = p2, height = 4, width = 6, dpi = 300)
saveRDS(so_ff_7000rpc, file.path(output_dir, "so_ff_7000rpc.rds"))


# --- Target 3: ON Dataset ---

cds_on <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-2-TRTL-ON_processed_cds_object.rds"))
# Note: this already has the cell type annotation if you just want to run this code and verify annotations, but use the output of fig1_7krpc for the preannotation cds

cds_on <- cds_on[,colData(cds_on)$n.umi > 70]
cds_on <- detect_genes(cds_on)

exprs_cds_on <- exprs(cds_on)
metadata_cds_on <- colData(cds_on)
so_on <- CreateSeuratObject(exprs_cds_on, meta.data = metadata_cds_on %>% as.data.frame(), assay = "RNA")
so_on <- NormalizeData(so_on)
so_on <- FindVariableFeatures(so_on, selection.method = "vst", nfeatures = 700)
so_on <- ScaleData(so_on)
so_on <- RunPCA(so_on, npcs = 20, verbose = F)
so_on <- RunUMAP(so_on, reduction = "pca", dims = 1:20, verbose = F)

so_on <- SCTransform(so_on)
so_on <- RunPCA(so_on, assay = "SCT")
so_on <- RunUMAP(so_on, dims = 1:15)

so.anchors <- FindTransferAnchors(
    reference = so_dt,
    query = so_on,
    dims = 1:15,
    reference.assay = "RNA", 
    query.assay = "RNA",
    features = unique(c(probe_genes, VariableFeatures(so_on, assay = "SCT")[1:200])),
    reference.reduction = "integrated.dr" 
)
predictions <- TransferData(anchorset = so.anchors, refdata = so_dt$major_cell_types, dims = 1:15)
so_on <- AddMetaData(so_on, metadata = predictions)

p3 <- DimPlot(so_on, group.by = "predicted.id") + facet_wrap("predicted.id")
ggsave(file.path(output_dir, "UMAPs/so_on_predicted.png"), plot = p3, height = 4, width = 6, dpi = 300)
saveRDS(so_on, file.path(output_dir, "so_on_final_22_dim.rds"))


# --- Target 4: dT 700 RPC ---
cds_dt_700rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_700rpc_cds_object.rds"))

# Note: this already has the cell type annotation if you just want to run this code and verify annotations, but use the output of fig1_700rpc for the preannotation cds

exprs_cds_dt_700rpc <- exprs(cds_dt_700rpc)
metadata_cds_dt_700rpc <- colData(cds_dt_700rpc)
so_dt_700rpc <- CreateSeuratObject(exprs_cds_dt_700rpc, meta.data = metadata_cds_dt_700rpc %>% as.data.frame(), assay = "RNA")
so_dt_700rpc <- NormalizeData(so_dt_700rpc)
so_dt_700rpc <- FindVariableFeatures(so_dt_700rpc, selection.method = "vst", nfeatures = 700)
so_dt_700rpc[["RNA"]] <- split(so_dt_700rpc[["RNA"]], f = so_dt_700rpc$experiment)
so_dt_700rpc <- SCTransform(so_dt_700rpc, variable.features.n = 500)
so_dt_700rpc <- RunPCA(so_dt_700rpc, assay = "SCT")
so_dt_700rpc <- IntegrateLayers(object = so_dt_700rpc, method = JointPCAIntegration, normalization.method = "SCT", verbose = F)
so_dt_700rpc <- FindNeighbors(so_dt_700rpc, reduction = "integrated.dr", dims = 1:20)
so_dt_700rpc <- RunUMAP(so_dt_700rpc, dims = 1:20, reduction = "integrated.dr")

so.anchors <- FindTransferAnchors(
    reference = so_dt,
    query = so_dt_700rpc,
    dims = 1:15,
    reference.assay = "RNA", 
    query.assay = "RNA",
    features = unique(c(probe_genes, VariableFeatures(so_dt, assay = "SCT")[1:500])),
    reference.reduction = "integrated.dr" 
)
predictions <- TransferData(anchorset = so.anchors, refdata = so_dt$major_cell_types, dims = 1:15)
so_dt_700rpc <- AddMetaData(so_dt_700rpc, metadata = predictions)
saveRDS(so_dt_700rpc, file.path(output_dir, "so_dt_700rpc_integrated_experiment.rds"))

# --- Target 5: TRTL 700 RPC (FF) ---
cds_ff_700rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_700rpc_cds_object.rds"))

# Note: this already has the cell type annotation if you just want to run this code and verify annotations, but use the output of fig1_700rpc for the preannotation cds

exprs_cds_ff_700rpc <- exprs(cds_ff_700rpc)
metadata_cds_ff_700rpc <- colData(cds_ff_700rpc)
so_ff_700rpc <- CreateSeuratObject(exprs_cds_ff_700rpc, meta.data = metadata_cds_ff_700rpc %>% as.data.frame(), assay = "RNA")
so_ff_700rpc <- NormalizeData(so_ff_700rpc)
so_ff_700rpc <- FindVariableFeatures(so_ff_700rpc, selection.method = "vst", nfeatures = 700)
so_ff_700rpc[["RNA"]] <- split(so_ff_700rpc[["RNA"]], f = so_ff_700rpc$experiment)
so_ff_700rpc <- SCTransform(so_ff_700rpc, variable.features.n = 700)
so_ff_700rpc <- RunPCA(so_ff_700rpc, assay = "SCT")
so_ff_700rpc <- IntegrateLayers(object = so_ff_700rpc, method = JointPCAIntegration, normalization.method = "SCT", verbose = F)
so_ff_700rpc <- FindNeighbors(so_ff_700rpc, reduction = "integrated.dr", dims = 1:15)
so_ff_700rpc <- RunUMAP(so_ff_700rpc, dims = 1:15, reduction = "integrated.dr")

so.anchors <- FindTransferAnchors(
    reference = so_dt,
    query = so_ff_700rpc,
    dims = 1:15,
    reference.assay = "RNA", 
    query.assay = "RNA",
    features = unique(c(probe_genes, VariableFeatures(so_dt, assay = "SCT")[1:500])), 
    reference.reduction = "integrated.dr"
)
predictions <- TransferData(anchorset = so.anchors, refdata = so_dt$major_cell_types, dims = 1:15)
so_ff_700rpc <- AddMetaData(so_ff_700rpc, metadata = predictions)
saveRDS(so_ff_700rpc, file.path(output_dir, "so_ff_700rpc.rds"))

# Update CDS with predictions
colData(cds_dt_7000rpc)$predicted.id <- so_dt_7000rpc$predicted.id
colData(cds_ff_7000rpc)$predicted.id <- so_ff_7000rpc$predicted.id
colData(cds_on)$predicted.id <- so_on$predicted.id
colData(cds_dt_700rpc)$predicted.id <- so_dt_700rpc$predicted.id
colData(cds_ff_700rpc)$predicted.id <- so_ff_700rpc$predicted.id

save_monocle_objects(cds_dt_7000rpc, file.path(output_dir, "cds/cds_dt_7000rpc_cell_type_annotated"))
save_monocle_objects(cds_ff_7000rpc, file.path(output_dir, "cds/cds_ff_7000rpc_cell_type_annotated"))
save_monocle_objects(cds_on, file.path(output_dir, "cds/cds_on_cell_type_annotated"))
save_monocle_objects(cds_dt_700rpc, file.path(output_dir, "cds/cds_dt_700rpc_cell_type_annotated"))
save_monocle_objects(cds_ff_700rpc, file.path(output_dir, "cds/cds_ff_700rpc_cell_type_annotated"))
