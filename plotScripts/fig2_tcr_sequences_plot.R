# Analyze and visualize TCR (TRA/TRB) sequence detection summaries.

library(dplyr)
library(ggplot2)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Load TCR Summary Data
# --------------------------------------------------------------------------------
TCR_summary <- readRDS(file.path(input_dir, "TCR_summary.RDS"))

# Create a summary dataframe to count cells based on TRA/TRB presence
plot_data <- TCR_summary %>%
  mutate(
    TRA_present = !is.na(TRA),
    TRB_present = !is.na(TRB),
    category = case_when(
      TRA_present & TRB_present ~ "TRA & TRB",
      TRA_present ~ "TRA Only",
      TRB_present ~ "TRB Only"
    ) 
  ) %>%
  dplyr::filter(!is.na(category))%>%
  dplyr::group_by(category) %>%
  dplyr::summarise(n = n())  

# --------------------------------------------------------------------------------
# Plot Bar Chart: TCR Categories
# --------------------------------------------------------------------------------
ggplot(plot_data, aes(x = category, y = n, fill = category)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25) +
  labs(
    title = "TCRs across 7,197\ncells (9.2% of total)",
    y = "Number of Cells",
    fill = "Category") +
  scale_fill_manual(values = c("TRA & TRB" = "maroon",
                               "TRA Only" = "orange",
                               "TRB Only" = "skyblue2")) +
  #theme_minimal() +
  monocle3:::monocle_theme_opts() +
  theme(
    text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    #plot.title = element_text(size = 9, face = "bold"),
    legend.position = "none")

ggsave(file.path(output_dir, "Proportion_of_cells_TCRA_TCRB.png"),height = 1.75, width = 2,dpi=900)

# --------------------------------------------------------------------------------
# Plot TRB Frequencies (Pie/Polar)
# --------------------------------------------------------------------------------
trb_counts <- TCR_summary %>%
  filter(!is.na(TRB)) %>%
  dplyr::count(TRB, sort = TRUE)
top3_trb <- trb_counts %>% slice_max(n, n = 3, with_ties = FALSE) %>% pull(TRB)
trb_counts <- trb_counts %>%
  mutate(fill_var = ifelse(TRB %in% top3_trb, TRB, "Other"))
# Truncate names for legend only
legend_labels <- c(setNames(paste0(substr(top3_trb, 1, 10), "..."), top3_trb), "Other" = "Other")
trb_counts$fill_var <- factor(trb_counts$fill_var, levels = c(top3_trb, "Other"))

p_trb <- ggplot(trb_counts, aes(x = "", y = n, fill = fill_var)) +
  geom_col(color = "grey80") +
  geom_bar(stat = "identity", width = 1, color = NA)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3, "Dark2"), "grey80")) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.25,"line")) +
  labs(title = "TRB", fill = "Top TCRβ sequences")

print(p_trb)
ggsave(file.path(output_dir, "TRB_sequences.png"),height = 0.8, width = 4.1,dpi=900)

# --------------------------------------------------------------------------------
# Plot TRA Frequencies (Pie/Polar)
# --------------------------------------------------------------------------------
tra_counts <- TCR_summary %>%
  filter(!is.na(TRA)) %>%
  dplyr::count(TRA, sort = TRUE)
top3_tra <- tra_counts %>% slice_max(n, n = 3, with_ties = FALSE) %>% pull(TRA)
tra_counts <- tra_counts %>%
  mutate(fill_var = ifelse(TRA %in% top3_tra, TRA, "Other"))
legend_labels <- c(setNames(paste0(substr(top3_tra, 1, 10), "..."), top3_tra), "Other" = "Other")
tra_counts$fill_var <- factor(tra_counts$fill_var, levels = c(top3_tra, "Other"))

p_tra <- ggplot(tra_counts, aes(x = "", y = n, fill = fill_var)) +
  geom_col(color = "grey80") +
  geom_bar(stat = "identity", width = 1, color = NA)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3, "Dark2"), "grey80")) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.25,"line")) +
  labs(title = "TRA", fill = "Top TCRα sequences")

print(p_tra)
ggsave(file.path(output_dir, "TRA_sequences.png"),height = 0.8, width = 4.1,dpi=900)


