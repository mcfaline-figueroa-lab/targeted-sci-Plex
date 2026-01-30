suppressPackageStartupMessages({
  library(devtools)
  library(parallel)
  library(ggplot2)
  library(ggridges)
  library(gridExtra)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(piano)
  library(UpSetR)
  library(DelayedArray)
  library(Seurat)
  library(monocle3)
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

# Create subdirectories for outputs
QC_plots_dir <- file.path(output_dir, "QC_plots")
marker_plots_dir <- file.path(output_dir, "Marker_plots")
umap_supp_dir <- file.path(output_dir, "UMAPs/markers_supp")

dir.create(QC_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(marker_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_supp_dir, recursive = TRUE, showWarnings = FALSE)


# --- Load Data ---

cds_path <- file.path(input_dir, "sci-Plex-mouse-brain_cds_precell_prehash.RDS")
if (!file.exists(cds_path)) stop("Input file not found: ", cds_path)
cds <- readRDS(cds_path)
names(colData(cds))[names(colData(cds)) == "Cell"] <- "cell"


# --- Cell Metadata Processing ---

process_metadata <- function(cds) {
  cols <- colData(cds) %>% 
    as.data.frame() %>%
    mutate(cell_id = cell) %>%
    tidyr::separate(col = cell_id, into = c('P7', 'P5', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
    mutate(RT = paste('RT_BC_', n1, sep = '')) %>%
    mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
    dplyr::select(-n1, -n2)
  
  return(cols)
}

cols <- process_metadata(cds)

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

cds <- process_cds_with_cols(cols, cds)

# --- Filtering & QC ---

mouse_genes <- rowData(cds)[grepl("ENSMUSG",rowData(cds)$id),]$id
cds <- cds[mouse_genes,]
cds <- cds[,colData(cds)$n.umi >= 100]

# QC Plots
p_umi <- ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(hash_umis)) %>% dplyr::mutate(rank = row_number()),
       aes(x = log10(rank), y = log10(hash_umis))) +
  geom_point() +
  geom_hline(yintercept = log10(10)) +
  monocle3:::monocle_theme_opts()
ggsave(file.path(scale_plots_dir, "Hash_UMIs_per_cell.png"), plot = p_umi, dpi = 600, width = 5, height = 3)

cds <- cds[,colData(cds)$hash_umis >= 5]

p_ratio <- ggplot(colData(cds) %>% as.data.frame(),
       aes(x = log10(top_to_second_best_ratio))) +
  geom_density() +
  geom_vline(xintercept = log10(2.5)) +
  monocle3:::monocle_theme_opts()
ggsave(file.path(scale_plots_dir, "Top_to_second_best_ratio.png"), plot = p_ratio, dpi = 600, width = 3, height = 3)

cds <- cds[,colData(cds)$top_to_second_best_ratio >= 2]

# --- Preprocessing ---

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)
expressed_genes <- rowData(cds) %>% as.data.frame() %>% filter(num_cells_expressed >= 50) %>% pull(id) %>% unique()
cds <- preprocess_cds(cds,
                      num_dim = 15,
                      use_genes = expressed_genes)

cds <- align_cds(cds, alignment_group = "oligo") # Note: alignment_group 'oligo' might need to exist in colData or be a placeholder from original script
cds <- reduce_dimension(cds,
                        preprocess_method = "Aligned",
                        umap.min_dist = 0.1,
                        umap.n_neighbors = 20)

colData(cds)$UMAP1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP2 <- reducedDims(cds)[["UMAP"]][,2]

cds <- cluster_cells(cds, reduction_method = "Aligned", resolution = 2e-3)
colData(cds)$Cluster <- monocle3::clusters(cds, reduction_method = "Aligned")

# Save Initial UMAP
p_umap <- plot_cells(cds, color_cells_by = "Cluster")
ggsave(file.path(output_dir, "initial_umap.png"), plot = p_umap)


# --- Marker Analysis ---

# Specific Genes Plot
p_genes <- plot_cells(cds, genes = c("Gfap","Olig2","Tubb3","Map2","Trem2","Tyrobp"), label_cell_groups = FALSE) +
  theme(legend.position = "right") +
  viridis::scale_color_viridis(option = "magma")
ggsave(file.path(output_dir, "selected_markers_umap.png"), plot = p_genes, width = 8, height = 6)


# Top Markers
marker_test_res <- top_markers(cds, group_cells_by="Cluster", 
                               reference_cells=1000, cores=1) # Reduced cores for safety

# QC Dist Plots for Markers
p_spec <- ggplot(marker_test_res, aes(x = specificity)) + 
  geom_density() + geom_vline(xintercept = 0.4) + monocle3:::monocle_theme_opts()
ggsave(file.path(scale_plots_dir, "marker_specificity.png"), plot = p_spec)

p_frac <- ggplot(marker_test_res, aes(x = fraction_expressing)) + 
  geom_density() + geom_vline(xintercept = 0.5) + monocle3:::monocle_theme_opts()
ggsave(file.path(scale_plots_dir, "marker_fraction.png"), plot = p_frac)

p_pseudo <- ggplot(marker_test_res, aes(x = pseudo_R2)) + 
  geom_density() + geom_vline(xintercept = 0.2) + monocle3:::monocle_theme_opts()
ggsave(file.path(scale_plots_dir, "marker_pseudoR2.png"), plot = p_pseudo)


top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10, specificity > 0.4, marker_test_q_value < 0.05) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

p_top <- plot_genes_by_group(cds,
                    c(top_specific_marker_ids),
                    group_cells_by="Cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3) +
  viridis::scale_color_viridis(option = "plasma") +
  theme(text =  element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(0.6,"line"), 
        legend.key.height = unit(0.6,"line")) +
  xlab("PCA cluster") +
  coord_flip()
ggsave(file.path(marker_plots_dir, "Top_markers.png"), plot = p_top, dpi = 900, width = 3.75, height = 1.5)


# --- Seurat Processing ---

cds_exprs<-exprs(cds)
so_cds <- CreateSeuratObject(cds_exprs, 
                             meta.data = colData(cds) %>% as.data.frame(), 
                             assay = "RNA")

so_cds <- NormalizeData(so_cds, normalization.method = "LogNormalize", scale.factor = 10000)
so_cds <- FindVariableFeatures(so_cds, selection.method = "vst", nfeatures = 2000)
so_cds <- ScaleData(so_cds)
so_cds <- RunPCA(so_cds, npcs = 30, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(so_cds[["pca"]], dims = 1:10, nfeatures = 5)
# DimHeatmap(so_cds, dims = 1:15, cells = 500, balanced = TRUE) # Optional, can be verbose

#Positive:  ENSMUSG00000002107, ENSMUSG00000054728, ENSMUSG00000060924, ENSMUSG00000062991, ENSMUSG00000036815 
#Negative:  ENSMUSG00000035202, ENSMUSG00000037625, ENSMUSG00000061808, ENSMUSG00000032554, ENSMUSG00000022548 
pc1_genes<-c("ENSMUSG00000002107", "ENSMUSG00000054728", "ENSMUSG00000060924", "ENSMUSG00000062991", "ENSMUSG00000036815","ENSMUSG00000035202", "ENSMUSG00000037625", "ENSMUSG00000061808", "ENSMUSG00000032554", "ENSMUSG00000022548") 
pc2_genes<-c("ENSMUSG00000054728", "ENSMUSG00000062991", "ENSMUSG00000060924", "ENSMUSG00000047495", "ENSMUSG00000002107","ENSMUSG00000041261", "ENSMUSG00000004630", "ENSMUSG00000028222", "ENSMUSG00000002930", "ENSMUSG00000042671")
pc3_genes<-c("ENSMUSG00000020428", "ENSMUSG00000060882", "ENSMUSG00000027965", "ENSMUSG00000055333", "ENSMUSG00000029088","ENSMUSG00000041261", "ENSMUSG00000028222", "ENSMUSG00000004630", "ENSMUSG00000030102", "ENSMUSG00000090223" )
pc4_genes<-c("ENSMUSG00000031425", "ENSMUSG00000037625", "ENSMUSG00000032554", "ENSMUSG00000015090", "ENSMUSG00000006782", "ENSMUSG00000017491", "ENSMUSG00000019990", "ENSMUSG00000054728", "ENSMUSG00000020599", "ENSMUSG00000060882") 
pc5_genes<-c("ENSMUSG00000036815", "ENSMUSG00000046159", "ENSMUSG00000062257", "ENSMUSG00000035357", "ENSMUSG00000059187", "ENSMUSG00000031425", "ENSMUSG00000037625", "ENSMUSG00000032554", "ENSMUSG00000015090", "ENSMUSG00000017491" )
pc6_genes<-c("ENSMUSG00000031425", "ENSMUSG00000053310", "ENSMUSG00000037625", "ENSMUSG00000035357", "ENSMUSG00000032554", "ENSMUSG00000034813", "ENSMUSG00000096914", "ENSMUSG00000052551", "ENSMUSG00000062209", "ENSMUSG00000035681")
pc7_genes<-c("ENSMUSG00000031425", "ENSMUSG00000037625", "ENSMUSG00000032554", "ENSMUSG00000032841", "ENSMUSG00000006782","ENSMUSG00000007097", "ENSMUSG00000079018", "ENSMUSG00000005360", "ENSMUSG00000024140", "ENSMUSG00000026728")
pc8_genes<-c("ENSMUSG00000020396", "ENSMUSG00000022054", "ENSMUSG00000022055", "ENSMUSG00000021948", "ENSMUSG00000026179","ENSMUSG00000052387", "ENSMUSG00000022708", "ENSMUSG00000005360", "ENSMUSG00000052920", "ENSMUSG00000007097")
pc9_genes<-c("ENSMUSG00000022054", "ENSMUSG00000022055", "ENSMUSG00000020396", "ENSMUSG00000062209", "ENSMUSG00000051111", "ENSMUSG00000021948", "ENSMUSG00000059857", "ENSMUSG00000040624", "ENSMUSG00000004698", "ENSMUSG00000052551")
pc10_genes<-c("ENSMUSG00000022212", "ENSMUSG00000055540", "ENSMUSG00000020524", "ENSMUSG00000042581", "ENSMUSG00000022708", "ENSMUSG00000063626", "ENSMUSG00000035357", "ENSMUSG00000048915", "ENSMUSG00000062991", "ENSMUSG00000058571")

so_cds <- FindNeighbors(so_cds, dims = 1:10)
so_cds <- FindClusters(so_cds, resolution = 0.5)
# head(Idents(so_cds), 5)
so_cds_markers <- FindAllMarkers(so_cds, only.pos = TRUE)
so_cds_markers_2FC <- so_cds_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  pull(gene)

# marker_genes_expanded will be used if needed for subsetting later.
marker_genes <- c("Hprt1","Pgk1","B2m","Plp1","Mbp","Cldn11","Trf","Cnp","Olig1","Olig2","Sox10","Gfap","Aqp4","Sox9","Slc7a10", "Mfsd2a","Agt", "Atp1a2","Slco1c1","Gpc5","Aldh1l1","S100b","Slc1a3","Slc1a2","Tubb3","Map2","Mapt","Gad2","Slc17a7","Kcna1","Cntn1","Car10","Kcnd2","Pde7b","Gabra6","Kcnd2","Dlx2","Gab1","Arx","Slc32a1","Dlx6as-1","Pcp2","Car8","Calb1","Itpr1","Trem2","Tmem119","Aif1","Cd68","Cd45")
marker_genes_id <- rowData(marker_genes_sub)$id
marker_genes_expanded<-c(marker_genes_id,pc1_genes,pc2_genes,pc3_genes,pc4_genes,pc5_genes,pc6_genes,pc7_genes,pc8_genes,pc9_genes,pc10_genes,so_cds_markers_2FC)
marker_genes_sub <- cds[rowData(cds)$id %in% marker_genes_expanded,]
marker_genes_expanded <- rowData(marker_genes_sub)$id
marker_genes_sub <- preprocess_cds(marker_genes_sub,
                                   num_dim = 10,
                                   use_genes = rowData(marker_genes_sub)$id)


marker_genes_sub<-align_cds(marker_genes_sub, alignment_group = "oligo")
marker_genes_sub <- reduce_dimension(marker_genes_sub,
                                     preprocess_method = "Aligned",
                                     umap.min_dist = 0.1,
                                     umap.n_neighbors = 20)

# save the main CDS here
save_monocle_objects(cds = cds, file.path(output_dir, "cds"))

# --- Sub-analysis with Probe Genes ---

genes_143_path <- file.path(input_dir, "Filtered_143_gene_R_L_final.csv")
if(file.exists(genes_143_path)){
  genes_143<-read.csv(genes_143_path)
  genes_143$gene_short_name<- sapply(strsplit(genes_143$probe_id, "_"), function(x) x[2])
  probe_genes <- unique(genes_143$gene_id)
  
  # Reload CDS (pattern from original script) or just use current
  # cds <- load_monocle_objects(...) # We have it in memory
  
  # Combine probe genes with the PC genes and Seurat markers we derived earlier
  marker_genes_expanded <- unique(c(probe_genes, marker_genes_pc, so_cds_markers_2FC))
  
  marker_genes_sub <- cds[rowData(cds)$id %in% marker_genes_expanded,]
  marker_genes_sub <- preprocess_cds(marker_genes_sub, num_dim = 10, use_genes = rowData(marker_genes_sub)$id)
  marker_genes_sub <- align_cds(marker_genes_sub, alignment_group = "oligo")
  marker_genes_sub <- reduce_dimension(marker_genes_sub, preprocess_method = "Aligned", umap.min_dist = 0.1, umap.n_neighbors = 20)
  
  colData(marker_genes_sub)$Cluster <- colData(cds)$Cluster
  
  save_monocle_objects(cds = marker_genes_sub, file.path(output_dir, "cds_marker_genes_sub"))
  
  # Dot plots loop
  marker_list <- list(
    "Oligodendrocytes" = c("Plp1", "Mbp", "Olig2", "Cldn11"),
    "Astrocytes"       = c("Gfap", "Aqp4", "Aldh1l1", "S100b"),
    "General_Neurons"  = c("Tubb3", "Map2", "Mapt"),
    "Excitatory_Neu"   = c("Car10", "Kcnd2", "Pde7b", "Slc17a7"),
    "Granule_Neurons"  = c("Gabra6", "Kcnd2"),
    "Interneurons"     = c("Dlx2", "Slc32a1", "Sst", "Gad2"),
    "Purkinje_Cells"   = c("Pcp2", "Calb1", "Itpr1", "Car8"),
    "Microglia"        = c("Tmem119", "Aif1", "Trem2", "Cd68"),
    "Endothelial"      = c("Phactr1", "Nrg1")
  )
  
  for (set_name in names(marker_list)) {
    valid_genes <- intersect(marker_list[[set_name]], rowData(cds)$gene_short_name)
    if(length(valid_genes) > 0) {
      p <- plot_genes_by_group(cds, 
                               markers = valid_genes, 
                               group_cells_by = "Cluster", 
                               ordering_type = "none",
                               reduction_method = "UMAP",
                               norm_method = "log",
                               max.size = 8) +
        scale_color_viridis_c(option = "magma") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10)) +
        labs(title = paste("Marker Group:", set_name), subtitle = "Clusters 1-8", x = "Marker Genes", y = "Cluster ID")
      
      print(p)
      ggsave(filename = file.path(marker_plots_dir, paste0("DotPlot_4inheight", set_name, ".png")), 
             plot = p, width = 4.5, height = 4, dpi = 900)
    }
  }
} else {
  warning("Filtered_143_gene_R_L_final.csv not found in input_dir. Skipping sub-analysis.")
}
