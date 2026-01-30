library(destiny)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(monocle3)
library(scran)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(densityClust)

# Arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide input_dir and output_dir as arguments.")
}
input_dir <- args[1]
output_dir <- args[2]

if (!dir.exists(file.path(output_dir, "plots"))) {
  dir.create(file.path(output_dir, "plots"), recursive = TRUE)
}

# Read CDS
cds_hash <- load_monocle_objects(file.path(input_dir, "targeted-sci-Plex-Jurkat_processed_cds_object.rds"))

# Probes and conversion to SingleCellExperiment
probe_genes <- c("CD4","CD69","IL2RA","PDCD1", "CD8A","IL2","TRAC","PTPRC","CD3G")
probe_genes_id <- rowData(cds_hash) %>% as.data.frame() %>% filter(gene_short_name %in% probe_genes) %>% pull(id) %>% unique()

jurkat_exprs <- exprs(cds_hash)
jurkat_meta <- as.data.frame(colData(cds_hash))
jurkat_row <- as.data.frame(rowData(cds_hash))
sce_jurkat <- SingleCellExperiment(
  assays = list(counts = jurkat_exprs),
  colData = jurkat_meta,
  rowData = jurkat_row
)

sce_jurkat <- logNormCounts(sce_jurkat)
top_hvg_genes <- modelGeneVar(sce_jurkat) %>%
  as.data.frame() %>%
  arrange(desc(bio)) %>%
  head(200) %>%
  rownames()
sce_jurkat <- sce_jurkat[unique(c(probe_genes_id, top_hvg_genes)),]

# Diffusion Map
reducedDim(sce_jurkat, "PCA_monocle") <- reducedDims(cds_hash)[["PCA"]]
dm <- DiffusionMap(reducedDim(sce_jurkat, "PCA_monocle"), verbose=TRUE)

dm_cells <- names(optimal_sigma(dm))
sce_jurkat_subset <- sce_jurkat[, dm_cells]
oligo <- sce_jurkat_subset$oligo
treatment <- sce_jurkat_subset$treatment
dose <- sce_jurkat_subset$dose

rm(dm_cells)

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  treatment = factor(sce_jurkat_subset$treatment, 
                                     levels = c("Media","DMSO","Dynabeads","ConcanavalinA","PHA","PMA","PMA+PHA")),
                  Dose= dose)

custom_colors <- c(
  "Media" = "lightgray",
  "DMSO" = "darkgray",
  "Dynabeads" = "#1B9E77",
  "ConcanavalinA" = "#D95F02",
  "PHA" = "#7570B3",
  "PMA" = "#E7298A",
  "PMA+PHA" = "#66A61E"
)

# DC1 vs DC2 Plot
ggplot(tmp, aes(x = DC1, y = DC2, colour = treatment)) +
  geom_point(alpha = 0.5) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    text = element_text(size = 12),
    legend.justification = c("right", "bottom"),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 7),
    legend.spacing.y = unit(-0.5, 'cm'),
    axis.text.x = element_text(angle = 45, hjust = 1) 
    
  ) +
  theme(axis.line = element_line(color = 'black')) +
  coord_cartesian(xlim = c(-.01, .05), ylim = c(-.025, .025))+
  scale_color_manual(values = custom_colors)

ggsave(file.path(output_dir, "plots/DC1vs2.png"), height = 3, width = 3, dpi = 900)

# Pseudotime Calculations (for Bar Plots)
sce_jurkat_subset$pseud_dm1 <- rank(eigenvectors(dm)[,1])
df_new <- data.frame(pseud_dm1 = sce_jurkat_subset$pseud_dm1, 
                     treatment = factor(sce_jurkat_subset$treatment, 
                                        levels = c("Media","DMSO","Dynabeads","ConcanavalinA","PHA","PMA","PMA+PHA")),
                     dose = sce_jurkat_subset$dose)

# Visualization: Treatment Composition by Trajectory Stage (Bar Plots)
my_breaks <- c(-Inf, 40000, Inf)
my_labels <- c("0-40k",  ">40k")

# Apply the cuts to create the cluster column
df_new$DC1_cluster <- cut(df_new$pseud_dm1, 
                          breaks = my_breaks, 
                          labels = my_labels,
                          include.lowest = TRUE)

plot_data <- df_new %>%
  group_by(treatment, DC1_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(treatment) %>%
  mutate(pct = count / sum(count))

ggplot(plot_data, aes(x = treatment, y = pct, fill = DC1_cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_text(data = subset(plot_data, DC1_cluster == ">40k"), 
            aes(label = scales::percent(pct, accuracy = 1)), 
            position = position_fill(vjust = 0.1), 
            color = "white", size = 2.5) +
  scale_fill_brewer(palette = "Set2") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Treatment Composition Activation Stages",
       x = "Treatment",
       y = "Proportion",
       fill = "Activation stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

ggsave(file.path(output_dir, "plots/histogram_stages.png"), height = 3, width = 4, dpi = 900)
