# Generate UMAP plots for various conditions (dT 7k, TRTL 7k, ON, dT full) 
# and compares them against Allen Brain Atlas annotations.

library(dplyr)
library(ggplot2)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory

# setwd(output_dir) # Optional: set working directory to output directory if needed

# --------------------------------------------------------------------------------
# Load dT 7k dataset and plot UMAPs for cell types (Neuron, Oligo, Astro, Endo)
# --------------------------------------------------------------------------------
## dT 16k Seurat transfer cell types to DT 7k reads:
cds_dt_7000rpc<-readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_7krpc_cds_object.rds"))

ggplot() +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(predicted.id == "Neuron"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(predicted.id == "Oligodendrocyte"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(predicted.id == "Astrocyte"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(predicted.id == "Endothelial"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  theme_void()

ggsave(file.path(output_dir, "UMAP_dt_7k_rpc_final.png"),height=0.75,width=0.75,dpi=900)

colData(cds_dt_7000rpc)$predicted.id <- factor(colData(cds_dt_7000rpc)$predicted.id, levels = c("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial"))

plot_cells(cds_dt_7000rpc[,!is.na(colData(cds_dt_7000rpc)$predicted.id)], color_cells_by = "predicted.id", label_cell_groups = FALSE) +
  theme(legend.position = "right",
        text = element_text(size = 7),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual("Cell type", 
                     values = c("Neuron" = "grey90", "Oligodendrocyte" = "navy", "Astrocyte" = "firebrick2", "Endothelial" = "deepskyblue2"),
                     labels = c("Neuron" = "Neu", "Oligodendrocyte" = "Olig", "Astrocyte" = "Ast", "Endothelial" = "End"))
ggsave(file.path(output_dir, "UMAP_dt_7k_rpc_final_for_legend.png"),height=1.5,width=1.5,dpi=900)

# --------------------------------------------------------------------------------
# Plot UMAPs using Allen Brain Atlas major cell type groups (dT 7k)
# --------------------------------------------------------------------------------
##dt ALLEN
ggplot() +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Glutamatergic"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "GABAergic"), aes(x = UMAP1, y =UMAP2), color = "grey70", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "OPC-Oligo"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Astro-Epen"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Vascular"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Immune"), aes(x = UMAP1, y =UMAP2), color = "brown", size = 0.2, stroke = 0) +
  theme_void()
ggsave(file.path(output_dir, "UMAP_dt_7k_rpc_final_Allen.png"),height=0.75,width=0.75,dpi=900)

##dt ALLEN for legend
colData(cds_dt_7000rpc)$major_group_cell_type <- factor(colData(cds_dt_7000rpc)$major_group_cell_type, levels = c("GABAergic", "Glutamatergic", "Astro-Epen", "OPC-Oligo", "Vascular", "Immune"))

plot_cells(cds_dt_7000rpc[,!is.na(colData(cds_dt_7000rpc)$major_group_cell_type)], color_cells_by = "major_group_cell_type", label_cell_groups = FALSE) +
  theme(legend.position = "right",
        text = element_text(size = 7),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual("Cell group", 
                     values = c("Glutamatergic" = "grey90", "GABAergic" = "grey70", "OPC-Oligo" = "navy", "Astro-Epen" = "firebrick2", "Vascular" = "deepskyblue2", "Immune" = "brown"),
                     labels = c("Glutamatergic" = "Glut", "GABAergic" = "GAPA", "OPC-Oligo" = "OPC-Olig", "Astro-Epen" = "Ast-Epen", "Vascular" = "Vascular", "Immune" = "Immune"))
ggsave(file.path(output_dir, "UMAP_dt_7k_rpc_final_Allen_for_legend.png"),height=1.5,width=1.5,dpi=900)

# --------------------------------------------------------------------------------
# Load TRTL 7k dataset, filter Microglia, and plot UMAPs
# --------------------------------------------------------------------------------
###TRTL
## dT 16k Seurat transfer cell types to TRTL 7k reads:
cds_TRTL_7000rpc<-readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_7krpc_cds_object.rds"))

cds_TRTL_mic_filt <- cds_TRTL_7000rpc[, colData(cds_TRTL_7000rpc)$predicted.id != "Microglia"]
cds_TRTL_mic_filt <- cds_TRTL_7000rpc[, colData(cds_TRTL_7000rpc)$predicted.id != "Microglia"]

colData(cds_TRTL_mic_filt)$UMAP1 <- reducedDims(cds_TRTL_mic_filt)[["UMAP"]][,1]
colData(cds_TRTL_mic_filt)$UMAP2 <- reducedDims(cds_TRTL_mic_filt)[["UMAP"]][,2]

ggplot() +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Neuron"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Oligodendrocyte"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Astrocyte"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Endothelial"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  theme_void()
ggsave(file.path(output_dir, "UMAP_TRTL_7k_rpc_final.png"),height=0.75,width=0.75,dpi=900)

#Allen plots for TRTL
colData(cds_TRTL_7000rpc)$UMAP1 <- reducedDims(cds_TRTL_7000rpc)[["UMAP"]][,1]
colData(cds_TRTL_7000rpc)$UMAP2 <- reducedDims(cds_TRTL_7000rpc)[["UMAP"]][,2]
ggplot() +
  geom_point(data = colData(cds_TRTL_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Glutamatergic"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "GABAergic"), aes(x = UMAP1, y =UMAP2), color = "grey70", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "OPC-Oligo"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Astro-Epen"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Vascular"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_7000rpc) %>% as.data.frame() %>% filter(major_group_cell_type == "Immune"), aes(x = UMAP1, y =UMAP2), color = "brown", size = 0.2, stroke = 0) +
  theme_void()
ggsave(file.path(output_dir, "UMAP_TRTL_7k_rpc_final_Allen.png"),height=0.75,width=0.75,dpi=900)


# --------------------------------------------------------------------------------
# Plot UMAPs for the ON condition (mapped from 16k dT)
# --------------------------------------------------------------------------------
#----------------------------------Supp ON--------------#
# Note: assuming cds_on is loaded from user environment or previously loaded. 
# If it needs loading: cds_on <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-2-TRTL-ON_processed_cds_object.rds"))

cds_on_mic_filt <- cds_on[, colData(cds_on)$predicted.id != "Microglia"]
cds_on_mic_filt <- cds_on[, colData(cds_on)$predicted.id != "Microglia"]

colData(cds_on_mic_filt)$UMAP1 <- reducedDims(cds_on_mic_filt)[["UMAP"]][,1]
colData(cds_on_mic_filt)$UMAP2 <- reducedDims(cds_on_mic_filt)[["UMAP"]][,2]

ggplot() +
  geom_point(data = colData(cds_on_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Neuron"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_on_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Oligodendrocyte"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_on_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Astrocyte"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_on_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Endothelial"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  theme_void()

ggsave(file.path(output_dir, "Supp_ON_UMAP_final.png"),height=0.75,width=0.75,dpi=900)

# --------------------------------------------------------------------------------
# Plot UMAPs for the full dT dataset and individual experiments split from it
# --------------------------------------------------------------------------------
#----------------------------------Supp dT full reads--------------#
colData(cds_dt)$UMAP1 <- reducedDims(cds_dt)[["UMAP"]][,1]
colData(cds_dt)$UMAP2 <- reducedDims(cds_dt)[["UMAP"]][,2]

ggplot() +
  geom_point(data = colData(cds_dt) %>% as.data.frame() %>% filter(major_cell_types == "Neuron"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt) %>% as.data.frame() %>% filter(major_cell_types == "Oligodendrocyte"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt) %>% as.data.frame() %>% filter(major_cell_types == "Astrocyte"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt) %>% as.data.frame() %>% filter(major_cell_types == "Endothelial"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  theme_void()

ggsave(file.path(output_dir, "UMAP_dt_full_rpc_final.png"),height=0.75,width=0.75,dpi=900)
#----------------------------------Supp dT full reads split--------------#
cds_list <- split(colnames(cds_dt), colData(cds_dt)$experiment) %>% 
  lapply(function(cells) cds_dt[, cells])

cds_list <- lapply(names(cds_list), function(exp_name) {
  cds_subset <- cds_list[[exp_name]]
  
  cds_subset <- preprocess_cds(cds_subset, method = "PCA")
  cds_subset <- reduce_dimension(cds_subset, 
                                 preprocess_method = "PCA",
                                 umap.min_dist = 0.1,
                                 umap.n_neighbors = 20,
                                 umap.metric = "cosine")
  
  colData(cds_subset)$UMAP1 <- reducedDims(cds_subset)[["UMAP"]][,1]
  colData(cds_subset)$UMAP2 <- reducedDims(cds_subset)[["UMAP"]][,2]
  
  p <- ggplot() +
    geom_point(data = colData(cds_subset) %>% as.data.frame() %>% filter(major_cell_types == "Neuron"), 
               aes(x = UMAP1, y = UMAP2), color = "grey90", size = 0.2, stroke = 0) +
    geom_point(data = colData(cds_subset) %>% as.data.frame() %>% filter(major_cell_types == "Oligodendrocyte"), 
               aes(x = UMAP1, y = UMAP2), color = "navy", size = 0.2, stroke = 0) +
    geom_point(data = colData(cds_subset) %>% as.data.frame() %>% filter(major_cell_types == "Astrocyte"), 
               aes(x = UMAP1, y = UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
    geom_point(data = colData(cds_subset) %>% as.data.frame() %>% filter(major_cell_types == "Endothelial"), 
               aes(x = UMAP1, y = UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
    theme_void()
  
  ggsave(file.path(output_dir, paste0("UMAP_dt_", exp_name, "_rpc_final.png")), 
         plot = p, height = 0.75, width = 0.75, dpi = 900)
  
  return(cds_subset)
})

# --------------------------------------------------------------------------------
# Generate dot plots for marker genes across cell types
# --------------------------------------------------------------------------------
#------------------------------------------------#

#dot-plot condensed markers dt-full#
marker_list <- list(
  "Oligodendrocytes"= c("Plp1", "Mbp", "Olig2", "Cldn11"),
  "Astrocytes"= c("Gfap", "Aqp4", "Aldh1l1", "S100b"),
  "General_Neurons"= c("Tubb3", "Map2", "Mapt"),
  "Excitatory_Neu"= c("Car10", "Kcnd2", "Pde7b", "Slc17a7"),
  "Granule_Neurons"= c("Gabra6", "Kcnd2"),
  "Interneurons"= c("Dlx2", "Slc32a1", "Sst", "Gad2"),
  "Purkinje_Cells"= c("Pcp2", "Calb1", "Itpr1", "Car8"),
  "Microglia"= c("Tmem119", "Aif1", "Trem2", "Cd68"),
  "Endothelial" = c("Phactr1","Ngr1")
  
)

#------------------------------Marker plots------------------------#
# setwd("~/Documents/runs/2307_sci/Marker_plots/") # Commented out as we use output_dir
for (set_name in names(marker_list)) {
  
  # Filter genes present in this specific dataset
  valid_genes <- intersect(marker_list[[set_name]], rowData(cds)$gene_short_name)
  
  if(length(valid_genes) > 0) {
    p <- plot_genes_by_group(cds, 
                             markers = valid_genes, 
                             group_cells_by = "Cluster", 
                             ordering_type = "none",   # Key for keeping your order
                             reduction_method = "UMAP",
                             norm_method = "log",
                             max.size = 8) +
      scale_color_viridis_c(option = "magma") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10)) +
      labs(title = paste("Marker Group:", set_name),
           subtitle = "Clusters 1-8",
           x = "Marker Genes",
           y = "Cluster ID")
    
    print(p)
    ggsave(filename = file.path(output_dir, paste0("DotPlot_4inheight", set_name, ".png")), 
           plot = p, 
           device = "png", 
           width = 4.5, 
           height = 4, 
           units = "in", 
           dpi = 900)
  }
}


