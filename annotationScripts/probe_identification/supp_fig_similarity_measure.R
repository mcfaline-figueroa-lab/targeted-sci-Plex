#from objects from fig1_sciPlex_mouse_brain_gene_panel_identification.R
marker_genes_sub<-load_monocle_objects("marker_gene_sub.rds")
cds<-load_monocle_objects("~cds.rds")

marker_genes_sub<-cluster_cells(marker_genes_sub, reduction_method = "Aligned", resolution = 1.5e-3)
colData(marker_genes_sub)$Cluster_new <- monocle3::clusters(marker_genes_sub, reduction_method = "Aligned")
length(unique(colData(marker_genes_sub)$Cluster_new))
cluster_comparison <- table(
  Original = colData(marker_genes_sub)$Cluster, 
  New = colData(marker_genes_sub)$Cluster_new
)

cluster_prop <- prop.table(cluster_comparison, margin = 1)

library(pheatmap)
png("cluster_retention_heatmap_900dpi.png",
    width = 4, height = 3, units = "in", res = 900)
pheatmap(cluster_prop, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "ghostwhite", "red"))(100),
         main = "Cluster retention",
         xlab = "Marker subset clusters",
         ylab = "Full transcriptome clusters")
dev.off()
plot_cells(marker_genes_sub, color_cells_by = "Cluster",group_label_size = 4) +
  scale_color_brewer(palette = "Paired") +
  theme_void()
ggsave("marker_genes_supp_clusters.png",dpi=900,height=2.5,width=3)
plot_cells(marker_genes_sub, color_cells_by = "Cluster_new",group_label_size = 4) +
  scale_color_brewer(palette = "Paired") +
  theme_void()
ggsave("marker_genes_supp_clusters_new.png",dpi=900,height=2.5,width=3)

############### genes UMAP
marker_list <- list(
  "Oligodendrocytes" = c("Plp1", "Mbp", "Olig2", "Cldn11"),
  "Astrocytes"       = c("Gfap", "Aqp4", "Aldh1l1", "S100b"),
  "General_Neurons"  = c("Tubb3", "Map2", "Mapt"),
  "Excitatory_Neu"   = c("Car10", "Kcnd2", "Pde7b", "Slc17a7"),
  "Granule_Neurons"  = c("Gabra6", "Kcnd2"),
  "Interneurons"     = c("Dlx2", "Slc32a1", "Sst", "Gad2"),
  "Purkinje_Cells"   = c("Pcp2", "Calb1", "Itpr1", "Car8"),
  "Microglia"        = c("Tmem119", "Aif1", "Trem2", "Cd68")
)
marker_list <- list(
  "Endothelial" = c("Phactr1", "Nrg1")
)
#dir.create("UMAPs/markers_supp")
setwd("UMAPs/markers_supp")
for (cell_type in names(marker_list)) {
  
  # Filter genes present in this specific dataset
  valid_genes <- intersect(marker_list[[cell_type]], rowData(marker_genes_sub)$gene_short_name)
  
  if(length(valid_genes) > 0) {
    p <- plot_cells(marker_genes_sub, 
                    genes = valid_genes,
                    label_cell_groups = FALSE,
                    show_trajectory_graph = FALSE) +
      # Force horizontal layout for this specific group
      facet_grid(. ~ feature_label) + 
      scale_color_viridis_c(option = "magma") +
      theme_void() +
      theme(strip.text.x = element_text(size = 10, face = "bold"),
            panel.spacing = unit(1, "lines")) +
      labs(title = paste("Markers for:", cell_type))
    
    print(p)
    
    # Calculate width based on number of genes to keep UMAPs square
    # 4 inches per gene is a good rule of thumb
    dynamic_width <- length(valid_genes) * 4
    
    ggsave(filename = paste0("Marker_umaps_", cell_type, ".png"), 
           plot = p, 
           width = dynamic_width, 
           height = 5, 
           units = "in", 
           dpi = 900)
  }
}

data.frame(n_umi = colData(cds)$n.umi) %>%
  arrange(desc(n_umi)) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = n_umi)) +
  geom_line() +
  scale_x_log10() + 
  scale_y_log10() +
  theme_minimal() +
  labs(x = "Rank", y = "n.UMI", title = "Knee Plot")

