# Generate supplementary UMAP plots for 700rpc datasets (dT and TRTL).

library(dplyr)
library(ggplot2)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Plot UMAPs for dT 700rpc
# --------------------------------------------------------------------------------
## dT 16k Seurat transfer cell types to DT 700 reads:
cds_dt_700rpc<-readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_700rpc_cds_object.rds"))

ggplot() +
  geom_point(data = colData(cds_dt_700rpc) %>% as.data.frame() %>% filter(predicted.id == "Neuron"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_700rpc) %>% as.data.frame() %>% filter(predicted.id == "Oligodendrocyte"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_700rpc) %>% as.data.frame() %>% filter(predicted.id == "Astrocyte"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_dt_700rpc) %>% as.data.frame() %>% filter(predicted.id == "Endothelial"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  theme_void()

ggsave(file.path(output_dir, "UMAP_dt_700_rpc_final.png"),height=0.75,width=0.75,dpi=900)

colData(cds_dt_700rpc)$predicted.id <- factor(colData(cds_dt_700rpc)$predicted.id, levels = c("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial"))

plot_cells(cds_dt_700rpc[,!is.na(colData(cds_dt_700rpc)$predicted.id)], color_cells_by = "predicted.id", label_cell_groups = FALSE) +
  theme(legend.position = "right",
        text = element_text(size = 7),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual("Cell type", 
                     values = c("Neuron" = "grey90", "Oligodendrocyte" = "navy", "Astrocyte" = "firebrick2", "Endothelial" = "deepskyblue2"),
                     labels = c("Neuron" = "Neu", "Oligodendrocyte" = "Olig", "Astrocyte" = "Ast", "Endothelial" = "End"))
ggsave(file.path(output_dir, "UMAP_dt_700_rpc_final_for_legend.png"),height=1.5,width=1.5,dpi=900)


plot_cells(cds_dt_700rpc[,!is.na(colData(cds_dt_700rpc)$major_group_cell_type)], color_cells_by = "major_group_cell_type", label_cell_groups = FALSE) +
  theme(legend.position = "right",
        text = element_text(size = 7),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line")) +
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual("Cell group", 
                     values = c("Glutamatergic" = "grey90", "GABAergic" = "grey70", "OPC-Oligo" = "navy", "Astro-Epen" = "firebrick2", "Vascular" = "deepskyblue2", "Immune" = "brown"),
                     labels = c("Glutamatergic" = "Glut", "GABAergic" = "GAPA", "OPC-Oligo" = "OPC-Olig", "Astro-Epen" = "Ast-Epen", "Vascular" = "Vascular", "Immune" = "Immune"))
ggsave(file.path(output_dir, "UMAP_dt_700_rpc_final_Allen_for_legend.png"),height=1.5,width=1.5,dpi=900)

# --------------------------------------------------------------------------------
# Plot UMAPs for TRTL 700rpc (TRTL 1k Reads)
# --------------------------------------------------------------------------------
###TRTL
## dT 16k Seurat transfer cell types to TRTL 700 reads:
cds_TRTL_700rpc<-readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_700rpc_cds_object.rds"))

cds_TRTL_mic_filt <- cds_TRTL_700rpc[, colData(cds_TRTL_700rpc)$predicted.id != "Microglia"]
cds_TRTL_mic_filt <- cds_TRTL_700rpc[, colData(cds_TRTL_700rpc)$predicted.id != "Microglia"]

colData(cds_TRTL_mic_filt)$UMAP1 <- reducedDims(cds_TRTL_mic_filt)[["UMAP"]][,1]
colData(cds_TRTL_mic_filt)$UMAP2 <- reducedDims(cds_TRTL_mic_filt)[["UMAP"]][,2]

ggplot() +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Neuron"), aes(x = UMAP1, y =UMAP2), color = "grey90", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Oligodendrocyte"), aes(x = UMAP1, y =UMAP2), color = "navy", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Astrocyte"), aes(x = UMAP1, y =UMAP2), color = "firebrick2", size = 0.2, stroke = 0) +
  geom_point(data = colData(cds_TRTL_mic_filt) %>% as.data.frame() %>% filter(predicted.id == "Endothelial"), aes(x = UMAP1, y =UMAP2), color = "deepskyblue2", size = 0.2, stroke = 0) +
  theme_void()
ggsave(file.path(output_dir, "UMAP_TRTL_700_rpc_final.png"),height=0.75,width=0.75,dpi=900)


