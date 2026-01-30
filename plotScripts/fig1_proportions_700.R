# Compare cell type proportions for 700rpc datasets (dT, TRTL) and relevant controls.

library(dplyr)
library(ggplot2)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

proportion_data <- readRDS(file.path(input_dir, "proportion_data_700rpc.RDS")) # Renamed from 1000.RDS
proportion_data$Dataset <- sapply(proportion_data$dataset, function(x){
  if(x == "DT:16k RPC")return("dT: 16k rpc")
  if(x == "DT:1k RPC")return("dT: 700 rpc")
  if(x == "FF:1k RPC")return("TRTL: 700 rpc")
})

custom_palette_fig1 <- (c("dT: 16k rpc"="orange",
                          "dT: 700 rpc" = "navyblue",
                          "TRTL: 700 rpc" = "skyblue2"))

# --------------------------------------------------------------------------------
# Plot Cell Type Proportions
# --------------------------------------------------------------------------------
ggplot(proportion_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell type", y = "Proportion of cells", fill = "Dataset") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.5,"line")) 
ggsave(file.path(output_dir, "Cell_type_proportions_full_vs_downsampled_700_rpc.png"),dpi=900,height = 1.75,width=2.5)

# --------------------------------------------------------------------------------
# Plot Astrocytes Proportions (Main + Inset)
# --------------------------------------------------------------------------------
astrocyte_data <- proportion_data %>% filter(cell_type == "Astrocyte")
p_astro_proportions <- ggplot(astrocyte_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell type", y = "Proportion of cells", fill = "Dataset") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.5,"line")) 
print(p_astro_proportions)
ggsave(file.path(output_dir, "Astro_proportions_full_vs_downsampled_700_rpc.png"),dpi=900,height = 1.75,width=2.5)
p_astro_proportions <- ggplot(astrocyte_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell type", y = "Proportion of cells", fill = "Dataset") +
scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(
    text = element_text(size = 7, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.5,"line")) 
ggsave(file.path(output_dir, "Astro_insert_proportions_full_vs_downsampled_700_rpc.png"),dpi=900,height = .9,width=.6)

# --------------------------------------------------------------------------------
# Plot Oligodendrocyte Proportions
# --------------------------------------------------------------------------------
oligo_data <- proportion_data %>% filter(cell_type == "Oligodendrocyte")
p_oligo_proportions <- ggplot(oligo_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell type", y = "Proportion of cells", fill = "Dataset") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(
    text = element_text(size = 9, color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.5,"line")) 
print(p_oligo_proportions)
ggsave(file.path(output_dir, "Oligo_proportions_full_vs_downsampled_700_rpc.png"),dpi=900,height = 1.75,width=2.5)


