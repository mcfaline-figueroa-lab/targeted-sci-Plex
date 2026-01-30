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
library(UpSetR)
library(DelayedArray)
library(monocle3)
library(readr)

# Reading in arguments from terminal
args = commandArgs(trailingOnly = T)
if (length(args) < 2) {
  stop("Please provide input_dir and output_dir as arguments.")
}
input_dir <- args[1]
output_dir <- args[2]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cds.pre.path = file.path(input_dir, "targeted-sci-Plex-Jurkat_cds_precell_prehash.RDS")
cds.pre <- readRDS(cds.pre.path)
names(colData(cds.pre))[names(colData(cds.pre)) == "Cell"] <- "cell"
RT.barcodes.path = file.path(input_dir, "sci-RNA-seq3_RT_shortdT_plate2.txt")
hash.path = file.path(input_dir, "targeted-sci-Plex-Jurkat_hashTable.out")

RT.barcodes <- read_table(file = RT.barcodes.path, col_names=c("RT", "Barcode"))
cols <- colData(cds.pre) %>% 
  as.data.frame() %>%
  mutate(cell_id = cell) %>%
  tidyr::separate(col = cell_id, into = c('P7', 'P5', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
  mutate(RT = paste('RT_', n1, sep = '')) %>%
  mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
  dplyr::select(-n1,-n2)

# Assigning hash to colData
hashTable <- read.table(hash.path, 
                        col.names = c("sample", "cell", "oligo", "hash_umis"))

hash_cells <- 1:length(hashTable$hash_umis)
hash_umis <- sort((hashTable$hash_umis),decreasing = TRUE)
high_umi_cells <- row.names(colData(cds.pre)[colData(cds.pre)$n.umi >= 100,])
hashTable %>% 
  filter(cell %in% high_umi_cells) %>% 
  dplyr::select(oligo) %>% 
  table()
png(file.path(output_dir, "hashknee.png"))
plot(log10(hash_cells),log10(hash_umis))
dev.off()

hashTable_summary <- hashTable %>%
  dplyr::group_by(cell) %>%
  dplyr::mutate(proportion = hash_umis/sum(hash_umis),total_hash_umis_per_cell = sum(hash_umis)) %>%
  dplyr::arrange(desc(proportion)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::filter(rank %in% c(1,2)) %>%
  dplyr::mutate(top_to_second_best_ratio = ifelse(sum(rank) > 1, 
                                                  proportion[rank == 1]/proportion[rank == 2], 1)) %>%
  dplyr::filter(rank == 1)

hashTable_summary$treatment <- sapply(hashTable_summary$oligo,function(x){stringr::str_split(x,pattern = "_")[[1]][1]})
hashTable_summary$dose <- sapply(hashTable_summary$oligo,function(x){stringr::str_split(x,pattern = "_")[[1]][2]})

cols <- cols %>% 
  left_join(hashTable_summary, by = c("cell" = "cell", "sample" = "sample")) %>%
  select(-Size_Factor, -proportion, -rank)

rownames(cols) <- cols$cell
cols <- DataFrame(cols)

identical(row.names(cols), row.names(colData(cds.pre)))
# assiging our cols dataframe as colData(cds)
colData(cds.pre) <- cols

# adding columns for log10(n.umi) and % mito
colData(cds.pre)$log10.umi <- colData(cds.pre)$n.umi %>% log10()

mt_genes <- rowData(cds.pre) %>% as.data.frame() %>%
  filter(grepl("mt-",gene_short_name) == TRUE)
mt_genes_id <- rownames(mt_genes)

colData(cds.pre)$percent_mito <- 100 * (colSums(exprs(cds.pre)[mt_genes_id,])/colSums(exprs(cds.pre)))

cells <- 1:length(colData(cds.pre)$cell)
umis <- sort(colData(cds.pre)$n.umi,decreasing = TRUE)
png(file.path(output_dir, "knee.png"))
plot(log10(cells),log10(umis))
dev.off()
cds.pre <- cds.pre[,colData(cds.pre)$n.umi>50]

cds.pre<-detect_genes(cds.pre)
expressed_genes_short_name <- rowData(cds.pre) %>% as.data.frame() %>% filter(num_cells_expressed >= .15*(length(colData(cds.pre)$cell))) %>% pull(gene_short_name) %>% unique()

expressed_genes <- rowData(cds.pre) %>% as.data.frame() %>% filter(num_cells_expressed >= .15*(length(colData(cds.pre)$cell))) %>% pull(id) %>% unique()
probe_genes<-c("CD4","CD69","IL2RA","PDCD1", "CD8A", "IL2","TRAC","PTPRC","CD3G")
probe_genes_id <- rowData(cds.pre) %>% as.data.frame() %>% filter(gene_short_name %in% probe_genes) %>% pull(id) %>% unique()
probe_genes_id <- rowData(cds_hash) %>% as.data.frame() %>% filter(gene_short_name %in% probe_genes) %>% pull(id) %>% unique()
length(expressed_genes)
cds.pre<-estimate_size_factors(cds.pre)
cds.pre <- preprocess_cds(cds.pre,
                          num_dim = 30,
                          use_genes = unique(c(expressed_genes,probe_genes_id)))

if (!dir.exists(file.path(output_dir, "QC_plots"))) {
  dir.create(file.path(output_dir, "QC_plots"))
}
plot_pc_variance_explained(cds.pre) +
  geom_vline(xintercept = 15)
ggsave(file.path(output_dir, "QC_plots/PC_variance_explained.png"), width = 2, height = 2, dpi = 900)

cds.pre <- preprocess_cds(cds.pre,
                          num_dim = 9,
                          use_genes = unique(c(expressed_genes,probe_genes_id)))
cds.pre <- reduce_dimension(cds.pre, 
                            max_components = 2,
                            preprocess_method = "PCA",
                            reduction_method = "UMAP",
                            umap.metric = "cosine",
                            umap.n_neighbors = 20,
                            umap.min_dist = 0.1,
                            umap.fast_sgd=FALSE,
                            cores=1,
                            verbose = T)
cds.pre<-estimate_size_factors(cds.pre)
plot_cells(cds.pre,scale_to_range = FALSE,graph_label_size = 20) 
cds.pre <- cluster_cells(cds.pre, resolution=5e-6,reduction_method="PCA")
colData(cds.pre)$Cluster <-clusters(cds.pre)
       length(unique(colData(cds.pre)$Cluster)) 
dim(cds.pre[,colData(cds.pre)$Cluster %in% "3"]) #45/78000 <1% of cells, outlier, filter out
cds.pre <- cds.pre[,!(colData(cds.pre)$Cluster %in% "3")]
dim(cds.pre)
plot_cells(cds_hash)+ facet_wrap(~treatment)+geom_density_2d(color = "black") 
if (!dir.exists(file.path(output_dir, "UMAP"))) {
  dir.create(file.path(output_dir, "UMAP"))
}
ggsave(file.path(output_dir, "UMAP/UMAP_bytreatment.png"),dpi=900,height=4,width=4)
cds_hash<-cds.pre[,!is.na(colData(cds.pre)$oligo)]
cds_hash<-cds_hash[,colData(cds_hash)$hash_umis>2]
cds_hash<-reduce_dimension(cds_hash, 
                 max_components = 2,
                 preprocess_method = "PCA",
                 reduction_method = "UMAP",
                 umap.metric = "cosine",
                 umap.n_neighbors = 20,
                 umap.min_dist = 0.05,
                 umap.fast_sgd=FALSE,
                 cores=1,
                 verbose = T)
#accounting for if there is only one hash (then top to second is one but it would be filtered out)
colData(cds_hash) <- DataFrame(as.data.frame(colData(cds_hash)) %>%
                                 mutate(top_to_second_best_ratio =
                                          if_else(hash_umis == total_hash_umis_per_cell,
                                                  10,
                                                  top_to_second_best_ratio)))

cds_hash<-cds_hash[,colData(cds_hash)$top_to_second_best_ratio>1.5]
plot_cells(cds_hash)+ facet_wrap(~oligo)+geom_density_2d(color = "black") 

ggsave(file.path(output_dir, "UMAP/by_oligo_clustered.png"),dpi=300)
plot_cells(cds_hash,genes = c(probe_genes), scale_to_range = FALSE, cell_size = .75) + 
  theme_minimal() +monocle3:::monocle_theme_opts() + 
  scale_color_viridis_c(option = "A") +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "none",
    #strip.text = element_blank(),  # Remove facet labels at the top
    plot.background = element_rect(fill = "white", color = NA))  # Set full plot background to white
    
ggsave(file.path(output_dir, "UMAP/by_oligo_probe_genes.png"),dpi=300)

plot_cells(cds.pre[,colData(cds.pre)$n.umi<2000],color_cells_by = "n.umi")+ facet_wrap(~P7)
cds.pre<-cds.pre[,colData(cds.pre)$n.umi<2000]
save_monocle_objects(cds.pre,file.path(output_dir, "cds_filt"))
cds_pre_hash<-load_monocle_objects(file.path(output_dir, "cds_filt"))
save_monocle_objects(cds_hash,file.path(output_dir, "cds_hash"))


#-------------------------combine TCR data

cols <- as.data.frame(colData(cds_hash)) %>%
  mutate(
    new_cell = paste(RT, Lig, P5, sep = "_") 
  )

TCR_summary<-readRDS(file.path(input_dir, "TCR_summary.RDS"))

dim(TCR_summary[TCR_summary$cell %in% cols$new_cell,]) #6751 
TCR_summary <- TCR_summary[TCR_summary$cell %in% cols$new_cell,]
cols_merge_TCR <- cols %>%
  left_join(TCR_summary, by = c("new_cell"="cell"))

rownames(cols_merge_TCR) <- cols_merge_TCR$cell
identical(rownames(colData(cds_hash)), row.names(cols_merge_TCR))

colData(cds_hash) <- DataFrame(cols_merge_TCR)
saveRDS(cds.pre, file.path(output_dir, "TCR_hash_merge_unfilt.cds"))
save_monocle_objects(cds_hash,file.path(output_dir, "cds_hash_with_TRA_TRB"))
cds.pre<-estimate_size_factors(cds.pre)

cds.pre <- preprocess_cds(cds.pre,
                          num_dim = 9,
                          use_genes = unique(c(expressed_genes,probe_genes_id)))
cds.pre <- reduce_dimension(cds.pre,
                            preprocess_method = "PCA",
                            umap.min_dist = .1,
                            umap.n_neighbors = 25)



Matrix::writeMM(t(exprs(cds_hash)), file = file.path(output_dir, "matrix.mtx"))
write.csv(as.data.frame(colData(cds_hash)), file = file.path(output_dir, "obs.csv"), row.names = TRUE)
write.csv(as.data.frame(rowData(cds_hash)), file = file.path(output_dir, "var.csv"), row.names = TRUE)


colData(cds_hash)$UMAP1 <- reducedDims(cds_hash)[["UMAP"]][,1]
colData(cds_hash)$UMAP2 <- reducedDims(cds_hash)[["UMAP"]][,2]
#save_monocle_objects(cds_hash,"~/Documents/runs/dt_combined_all_reads/cds_hash")


