# Compare UMI counts and marker gene expression between oligo-dT and TRTL methods.

library(dplyr)
library(ggplot2)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Load Data and Prepare Combined CDS (7k RPC)
# --------------------------------------------------------------------------------
# Load necessary objects (assuming they were saved previously or need to be loaded)
# Note: combine_cds needs cds objects. Loading them here if they aren't in environment.
cds_dt_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_7krpc_cds_object.rds"))
cds_TRTL_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_7krpc_cds_object.rds"))


cds_comb_7000<-combine_cds(list(cds_dt_7000rpc,cds_TRTL_7000rpc))
astrocyte_sub<-cds_comb_7000[,colData(cds_comb_7000)$predicted.id %in% "Astrocyte"]
colData(astrocyte_sub)$log.n.umi<-colData(astrocyte_sub)$log10.umi
#save_monocle_objects(astrocyte_sub,file.path(output_dir, "astrocyte_combined_7000rpc.RDS"))

colData(astrocyte_sub)$method <- sapply(colData(astrocyte_sub)$method, function(x){
  if(x == "dt")return("oligo-dT")
  if(x == "ff")return("TRTL")
})

# --------------------------------------------------------------------------------
# Plot UMI Counts for Astrocytes
# --------------------------------------------------------------------------------
colData(astrocyte_sub) %>%
  as.data.frame() %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = log10(n.umi), fill = method), size = 0.1, outlier.shape = NA) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values =  c("oligo-dT" = "dimgrey", "TRTL" = "deepskyblue2")) +
  theme(text = element_text(size = 9),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7)) +
  xlab("Method") +
  ylab("log10(UMIs)")
ggsave(file.path(output_dir, "Log10_astroscyte_UMIs.png"),dpi=900,height = 1.5,width= 0.8)

colData(astrocyte_sub) %>%
  as.data.frame() %>%
  group_by(method) %>%
  dplyr::summarize(med_umi = median(n.umi),
                   med_genes_expr = median(num_genes_expressed),
                   number_cells = n())
# method   med_umi med_genes_expr number_cells
# <chr>      <dbl>          <dbl>        <int>
# 1 TRTL         722            497         2298
# 2 oligo-dT    2167           1423         1235
#subsample to 720UMI for normalized size factor based on UMI
astrocyte_sub_exprs<-exprs(astrocyte_sub)
astrocyte_sub_exprs<-subsample_umis_to_target_median(astrocyte_sub_exprs, 720)
astrocyte_sub_UMI_subsampled<-new_cell_data_set(astrocyte_sub_exprs, as.data.frame(colData(astrocyte_sub)), as.data.frame(rowData(astrocyte_sub)))
astrocyte_sub_UMI_subsampled<-estimate_size_factors(astrocyte_sub_UMI_subsampled)
astrocyte_sub_UMI_subsampled<-detect_genes(astrocyte_sub_UMI_subsampled)
colData(astrocyte_sub_UMI_subsampled)$n.umi<-colSums(exprs(astrocyte_sub_UMI_subsampled))
colData(astrocyte_sub_UMI_subsampled) %>%
  as.data.frame() %>%
  group_by(method) %>%
  dplyr::summarize(med_umi = median(n.umi),
                   med_genes_expr = median(num_genes_expressed),
                   number_cells = n())

save_monocle_objects(astrocyte_sub_UMI_subsampled,file.path(output_dir, "astrocyte_sub_UMI_subsampled.rds")) # Updated path

#astrocyte_sub_UMI_subsampled <- readRDS(file.path(input_dir, "astrocyte_sub_UMI_subsampled.rds"))

colData(astrocyte_sub_UMI_subsampled)$method <- sapply(colData(astrocyte_sub_UMI_subsampled)$method, function(x){
  if(x == "dt")return("oligo-dT")
  if(x == "TRTL")return("TRTL")
  if(x == "ff")return("TRTL")
})

# --------------------------------------------------------------------------------
# Plot Expression of Astrocyte Markers
# --------------------------------------------------------------------------------
marker_genes_to_plot<-c("Gfap","Aqp4","Sox9","Slc7a10", "Mfsd2a", "Atp1a2","Slco1c1","Gpc5","Aldh1l1","S100b","Slc1a3","Slc1a2")

astrocyte_sub_UMI_subsampled_genes<- astrocyte_sub_UMI_subsampled[rowData(astrocyte_sub_UMI_subsampled)$gene_short_name %in% marker_genes_to_plot, ]
plot_percent_cells_positive(astrocyte_sub_UMI_subsampled_genes, 
                            group_cells_by = "method", 
                            min_expr = 1, 
                            normalize = TRUE, 
                            ncol = 6) + 
  ggplot2::scale_fill_manual("Method", values = c("oligo-dT" = "dimgrey", "TRTL" = "deepskyblue2")) +
  theme(
    text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "Astrocyte_markers_nopval.png"),height=1.75,width= 3.25,dpi=900)
plot_percent_cells_positive_pvals(astrocyte_sub_UMI_subsampled_genes, 
                            group_cells_by = "method", 
                            min_expr = 1, 
                            normalize = TRUE, 
                            ncol = 6) + 
  ggplot2::scale_fill_manual("Method", values = c("oligo-dT" = "dimgrey", "TRTL" = "deepskyblue2")) +
  theme(
    text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position = "none"
  )
ggsave(file.path(output_dir, "Astrocyte_markers.png"),height=1.75,width= 3.25,dpi=900)

# --------------------------------------------------------------------------------
# Oligodendrocyte Analysis (Supplementary Figure)
# --------------------------------------------------------------------------------
####Supp fig oligodendrocyte same method:
oligo_sub<-cds_comb_7000[,colData(cds_comb_7000)$predicted.id %in% "Oligodendrocyte"]
colData(oligo_sub)$log.n.umi<-colData(oligo_sub)$log10.umi
save_monocle_objects(oligo_sub,file.path(output_dir, "oligo_combined_7000rpc.RDS"))

# oligo_sub <- readRDS(file.path(input_dir, "oligo_sub.rds")) # Updated to input_dir. Note: ensure this file exists or is created.

colData(oligo_sub)$method <- sapply(colData(oligo_sub)$method, function(x){
  if(x == "dt")return("oligo-dT")
  if(x == "TRTL")return("TRTL")
  if(x == "ff")return("TRTL")
})

colData(oligo_sub) %>%
  as.data.frame() %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = log10(n.umi), fill = method), size = 0.1, outlier.shape = NA) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values =  c("oligo-dT" = "dimgrey", "TRTL" = "deepskyblue2")) +
  theme(text = element_text(size = 9),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7)) +
  xlab("Method") +
  ylab("log10(UMIs)")
ggsave(file.path(output_dir, "Log10_oligo_UMIs.png"),dpi=900,height = 1.5,width= 0.8)
colData(oligo_sub) %>%
  as.data.frame() %>%
  group_by(method) %>%
  dplyr::summarize(med_umi = median(n.umi),
                   med_genes_expr = median(num_genes_expressed),
                   number_cells = n())
# method   med_umi med_genes_expr number_cells
# <chr>      <dbl>          <int>        <int>
# 1 TRTL        1070            707         2039
# 2 oligo-dT    2343           1519         1099
oligo_sub_exprs<-exprs(oligo_sub)
oligo_sub_exprs<-subsample_umis_to_target_median(oligo_sub_exprs, 1000)
oligo_sub_UMI_subsampled<-new_cell_data_set(oligo_sub_exprs, as.data.frame(colData(oligo_sub)), as.data.frame(rowData(oligo_sub)))
oligo_sub_UMI_subsampled<-estimate_size_factors(oligo_sub_UMI_subsampled)
oligo_sub_UMI_subsampled<-detect_genes(oligo_sub_UMI_subsampled)
colData(oligo_sub_UMI_subsampled)$n.umi<-colSums(exprs(oligo_sub_UMI_subsampled))
colData(oligo_sub_UMI_subsampled) %>%
  as.data.frame() %>%
  group_by(method) %>%
  dplyr::summarize(med_umi = median(n.umi),
                   med_genes_expr = median(num_genes_expressed),
                   number_cells = n())

save_monocle_objects(oligo_sub_UMI_subsampled, file.path(output_dir, "oligo_sub_UMI_subsampled.rds"))

# oligo_sub_UMI_subsampled <- readRDS(file.path(input_dir, "oligo_sub_UMI_subsampled.rds"))

colData(oligo_sub_UMI_subsampled)$method <- sapply(colData(oligo_sub_UMI_subsampled)$method, function(x){
  if(x == "dt")return("oligo-dT")
  if(x == "TRTL")return("TRTL")
  if(x == "ff")return("TRTL")
})

# --------------------------------------------------------------------------------
# Plot Expression of Oligodendrocyte Markers
# --------------------------------------------------------------------------------
marker_genes_to_plot<-c("Plp1","Mbp","Cldn11","Trf","Cnp","Olig1","Olig2","Sox10")

oligo_sub_UMI_subsampled_genes<- oligo_sub_UMI_subsampled[rowData(oligo_sub_UMI_subsampled)$gene_short_name %in% marker_genes_to_plot, ]
plot_percent_cells_positive(oligo_sub_UMI_subsampled_genes, 
                            group_cells_by = "method", 
                            min_expr = 1, 
                            normalize = TRUE, 
                            ncol = 4) + 
  ggplot2::scale_fill_manual("Method", values = c("oligo-dT" = "dimgrey", "TRTL" = "deepskyblue2")) +
  theme(
    text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position = "none"
  )
ggsave(file.path(output_dir, "oligo_markers_no_pval.png"),height=1.75,width= 3.25,dpi=900)

plot_percent_cells_positive_pvals(oligo_sub_UMI_subsampled_genes, 
                                  group_cells_by = "method", 
                                  min_expr = 1, 
                                  normalize = TRUE, 
                                  ncol = 4) + 
  ggplot2::scale_fill_manual("Method", values = c("oligo-dT" = "dimgrey", "TRTL" = "deepskyblue2")) +
  theme(
    text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position = "none"
  )
ggsave(file.path(output_dir, "oligo_markers.png"),height=1.75,width= 3.25,dpi=900)
