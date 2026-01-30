# Analyze Jurkat activation trajectory and marker expression across pseudotime.

# library(destiny)
# library(Seurat)
# library(SingleCellExperiment)
library(scater)
# library(scuttle)
library(ggplot2)
library(dplyr)
library(monocle3)
# library(scran)
library(pheatmap)
library(RColorBrewer)
# library(densityClust)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Density Curve Analysis
# --------------------------------------------------------------------------------
df_plot_final <- readRDS(file.path(input_dir, "df_plot_final.RDS"))

# Map treatment names for plotting
df_plot_final$Treatment <- sapply(df_plot_final$Plot_Group, function(x){
  if(x == "Dynabeads") return("Dyna")
  if(x == "ConcanavalinA") return("ConA")
  if(x == "PHA") return("PHA")
  if(x == "PMA") return("PMA")
  if(x == "PMA+PHA") return("PMA+PHA")
  if(x == "Media Control") return("Media")
  if(x == "DMSO Control") return("DMSO")
  return(NA)
})

df_plot_final$Treatment <- factor(df_plot_final$Treatment , levels = c("Media", "DMSO", "Dyna", "ConA", "PHA", "PMA", "PMA+PHA"))

df_plot_final$Facet_Group <- sapply(df_plot_final$Facet_Group, function(x){
  if(x == "Dynabeads") return("Dyna")
  if(x == "ConcanavalinA") return("ConA")
  if(x == "PHA") return("PHA")
  if(x == "PMA") return("PMA")
  if(x == "PMA+PHA") return("PMA+PHA")
  return(NA)
})

custom_palette <- c(
  "Dyna" = "#1B9E77",
  "ConA" = "#D95F02",
  "PHA" = "#7570B3",
  "PMA" = "#E7298A",
  "PMA+PHA" = "#66A61E",
  "Media" = "gray60",
  "DMSO" = "gray80"
)

ggplot(df_plot_final, aes(x = pseud_dm1, color = Treatment, group = Treatment)) +
  geom_density(size = 0.5) +
  monocle3:::monocle_theme_opts() +
  theme(
    text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    legend.position = "none",
    legend.key.width = unit(0.2,"line"),
    legend.key.height = unit(0.25,"line")) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  scale_color_manual(values = custom_palette, name = "Treatment") +
  facet_grid(. ~ Facet_Group, scales = "free_y") +
  xlab("Activation Trajectory (pseudotime)") +
  ylab("Cell density")

ggsave(file.path(output_dir, "Jurkat_density_across_pseudotime_by_treatment.png"), dpi=900, height = 1.25, width = 3.5)

# --------------------------------------------------------------------------------
# Marker Expression Across Pseudotime
# --------------------------------------------------------------------------------
sce_jurkat_subset <- readRDS(file.path(input_dir, "sce_jurkat_subset.RDS"))

plotExpression(sce_jurkat_subset, c("IL2","IL2RA","PDCD1"), x = "pseud_dm1", 
               show_violin = FALSE, colour_by = "treatment",
               show_smooth = FALSE,ncol = 3) +
  scale_color_manual(values = custom_palette)+
  #theme_void() +
  monocle3:::monocle_theme_opts() +
  theme(
    text = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    legend.position = "none",
    strip.text = element_text(face = "italic")) +
  xlab("Activation Trajectory (pseudotime)") +
  ylab("Expression\n(logcounts)")

ggsave(file.path(output_dir, "Jurkat_marker_expression_across_pseudotime.png"), dpi=900,height = 1.25,width = 3.5)


