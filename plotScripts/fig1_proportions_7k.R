# Calculate and plot cell type proportions across different datasets 
# (dT 16k, dT 7k, TRTL 7k) and compares them with Allen Brain Atlas data.

library(dplyr)
library(ggplot2)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Load Data and Prepare Proportions (Seurat Label Transfer)
# --------------------------------------------------------------------------------
#####Fig 1 cell type proportions SEURAT, dT 16000 is reference (bottom is Allen)
# 
# Load necessary files (assuming they are needed to construct the dataframes)
cds_dt <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_16krpc_cds_object.rds"))
cds_dt_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_7krpc_cds_object.rds"))
cds_TRTL_7000rpc <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_TRTL_7krpc_cds_object.rds"))


dt_data <- as.data.frame(colData(cds_dt)) %>%
  mutate(dataset = "DT-16k-RPC") %>%
  rename(cell_type = major_cell_types) %>%
  select(cell_type, dataset)

dt_7000_data <- as.data.frame(colData(cds_dt_7000rpc)) %>%
  mutate(dataset = "DT-7k-RPC") %>%
  rename(cell_type = predicted.id) %>%
  select(cell_type, dataset)

TRTL_7000_data <- as.data.frame(colData(cds_TRTL_7000rpc)) %>%
  mutate(dataset = "TRTL-7k-RPC") %>%
  rename(cell_type = predicted.id) %>%
  select(cell_type, dataset)

list_of_metadata <- list(dt_data, dt_7000_data, TRTL_7000_data)
ordered_datasets <- c("DT-16k-RPC", "DT-7k-RPC", "TRTL-7k-RPC")

unified_data_long <- bind_rows(list_of_metadata) %>%
  mutate(
    cell_type = as.factor(cell_type),
    dataset = factor(dataset, levels = ordered_datasets)
  ) %>%
  select(dataset, cell_type) %>%
  drop_na(cell_type)

proportion_data <- unified_data_long %>%
  group_by(dataset) %>%
  mutate(Total_Cells = n()) %>%
  ungroup() %>%
  group_by(dataset, cell_type, Total_Cells) %>%
  summarise(Cell_Count = n(), .groups = 'drop') %>%
  mutate(Proportion = Cell_Count / Total_Cells) %>%
  ungroup()

ordered_types <- proportion_data %>%
  filter(dataset == "DT-16k-RPC") %>%
  arrange(desc(Proportion)) %>%
  pull(cell_type)
proportion_data$cell_type <- factor(proportion_data$cell_type, levels = ordered_types)

#proportion_data<-readRDS(file.path(input_dir, "proportion_data_7000.RDS"))

# Add a placeholder for microglia (not called in downsampled datasets) for plotting
proportion_data <- rbind(proportion_data %>% dplyr::select(-Total_Cells), 
                         data.frame(dataset = c("DT-7k-RPC","TRTL-7k-RPC"), 
                                    cell_type = c("Microglia","Microglia"),
                                    Cell_Count = c(0,0),
                                    Proportion = c(0,0)))

proportion_data$Dataset <- sapply(proportion_data$dataset, function(x){
  if(x == "DT-16k-RPC")return("dT:16k rpc")
  if(x == "DT-7k-RPC")return("dT:7k rpc")
  if(x == "TRTL-7k-RPC")return("TRTL:7k rpc")
})

custom_palette_fig1 <- (c("dT:16k rpc"="orange",
                          "dT:7k rpc" = "navyblue",
                          "TRTL:7k rpc" = "skyblue2"))

# --------------------------------------------------------------------------------
# Plot General Cell Type Proportions
# --------------------------------------------------------------------------------
ggplot(proportion_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Cell_type_proportions_full_vs_downsampled_Seurat_label_transfer.png"),
       dpi=900,height = 1.75,width=2.5)

# Plot Oligodendrocytes specifically
oligo_data <- proportion_data %>% filter(cell_type == "Oligodendrocyte")

ggplot(oligo_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Oligodendrocyte_proportions_sci_reference_for_supplement.png"),
       dpi=900,height = 1.75,width=2)

# --------------------------------------------------------------------------------
# Allen Brain Atlas Proportions
# --------------------------------------------------------------------------------
proportion_data_allen <- readRDS(file.path(input_dir, "proportion_data_7000_allen.RDS"))

proportion_data_allen$Dataset <- sapply(proportion_data_allen$dataset, function(x){
  if(x == "DT-16k-RPC")return("dT:16k rpc")
  if(x == "DT-7k-RPC")return("dT:7k rpc")
  if(x == "TRTL-7k-RPC")return("TRTL:7k rpc")
})

ggplot(proportion_data_allen, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Cell_type_proportions_full_vs_downsampled_Allen.png"),
       dpi=900,height = 1.75,width=2.5)

# Plot Astrocytes specifically (Allen)
astrocyte_data <- proportion_data_allen %>% filter(cell_type == "Astro-Epen")

ggplot(astrocyte_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Astrocyte_proportions_full_vs_downsampled_Allen.png"),
       dpi=900,height = 0.9,width=0.5)

ggplot(astrocyte_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Astrocyte_proportions_full_vs_downsampled_Allen_for_supplement.png"),
       dpi=900,height = 1.75,width=2)

# Plot Oligodendrocytes specifically (Allen)
oligodendrocyte_data <- proportion_data_allen %>% filter(cell_type == "OPC-Oligo")

ggplot(oligodendrocyte_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Oligodendrocyte_proportions_full_vs_downsampled_Allen_for_supplement.png"),
       dpi=900,height = 1.75,width=2)

# Plot Immune specifically (Allen)
mic_data <- proportion_data_allen %>% filter(cell_type == "Immune")
ggplot(mic_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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
ggsave(file.path(output_dir, "Immune_proportions_full_vs_downsampled_Allen_for_supplement.png"),
       dpi=900,height = 1.75,width=2)

# --------------------------------------------------------------------------------
# ON Condition Proportions
# --------------------------------------------------------------------------------

# Read files if not already loaded
# cds_dt <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-1-2-combined_processed_dT_16krpc_cds_object.rds"))
cds_on <- readRDS(file.path(input_dir, "sci-RNA-TRTL-mouse-2-TRTL-ON_processed_cds_object.rds"))

dt_data <- as.data.frame(colData(cds_dt)) %>%
  mutate(dataset = "DT-16k-RPC") %>%
  rename(cell_type = major_cell_types) %>%
  select(cell_type, dataset)

on_data <- as.data.frame(colData(cds_on)) %>%
  mutate(dataset = "TRTL-ON-500-RPC") %>%
  rename(cell_type = predicted.id) %>%
  select(cell_type, dataset)


list_of_metadata <- list(dt_data,on_data)
ordered_datasets <- c("DT-16k-RPC", "TRTL-ON-500-RPC")

unified_data_long <- bind_rows(list_of_metadata) %>%
  mutate(
    cell_type = as.factor(cell_type),
    dataset = factor(dataset, levels = ordered_datasets)
  ) %>%
  select(dataset, cell_type) %>%
  drop_na(cell_type)

proportion_data <- unified_data_long %>%
  group_by(dataset) %>%
  mutate(Total_Cells = n()) %>%
  ungroup() %>%
  group_by(dataset, cell_type, Total_Cells) %>%
  summarise(Cell_Count = n(), .groups = 'drop') %>%
  mutate(Proportion = Cell_Count / Total_Cells) %>%
  ungroup()

ordered_types <- proportion_data %>%
  filter(dataset == "DT-16k-RPC") %>%
  arrange(desc(Proportion)) %>%
  pull(cell_type)
proportion_data$cell_type <- factor(proportion_data$cell_type, levels = ordered_types)

#proportion_data<-readRDS("~/Desktop/Figures_targeted/proportion_data_7000.RDS")

# Add a placeholder for microglia (not called in downsampled datasets) for plotting
proportion_data <- rbind(proportion_data %>% dplyr::select(-Total_Cells), 
                         data.frame(dataset = c("TRTL-ON-500-RPC"), 
                                    cell_type = c("Microglia"),
                                    Cell_Count = c(0),
                                    Proportion = c(0)))

proportion_data$Dataset <- sapply(proportion_data$dataset, function(x){
  if(x == "DT-16k-RPC")return("dT:16k rpc")
  if(x == "TRTL-ON-500-RPC")return("TRTL-ON:500 rpc")
  
})

custom_palette_fig1 <- (c("dT:16k rpc"="orange",
                          "TRTL-ON:500 rpc" = "skyblue4"))

ggplot(proportion_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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
ggsave(file.path(output_dir, "Cell_type_proportions_full_vs_TRTL-ON.png"),
       dpi=900,height = 1.75,width=2.5)

oligo_data <- proportion_data %>% filter(cell_type == "Oligodendrocyte")

ggplot(oligo_data, aes(x = cell_type, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette_fig1) +
  labs(x = "Cell Type", y = "Proportion of cells", fill = "Dataset") +
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

ggsave(file.path(output_dir, "Oligodendrocyte_proportions_ON_sci_reference_for_supplement.png"),
       dpi=900,height = 1.75,width=2)

# --------------------------------------------------------------------------------
#Pearson
# --------------------------------------------------------------------------------

print(paste0("Pearson Correlation (DT-7k vs 16k): ", round(cor_dt_7k, 4)))
#"Pearson Correlation (DT-7k vs 16k): 0.9992"
print(paste0("Pearson Correlation (TRTL-7k vs 16k): ", round(cor_trtl_7k, 4)))
#"Pearson Correlation (TRTL-7k vs 16k): 0.9996"
