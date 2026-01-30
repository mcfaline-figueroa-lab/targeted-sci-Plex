library(dplyr)
library(ggplot2)
library(monocle3)

setwd("~/Dropbox (Personal)/McFaline-Figueroa Lab/Publications/targeted_sci-Plex/Figures/")

cds_hash<-readRDS("~/Dropbox (Personal)/McFaline-Figueroa Lab/Publications/targeted_sci-Plex/Figures/cds_object.rds")

plot_data <- colData(cds_hash) %>%
  as.data.frame() %>%
  mutate(
    has_TRA = !is.na(TRA) & (TRA != "") & (TRA != 0),
    has_TRB = !is.na(TRB) & (TRB != "") & (TRB != 0),
    TCR_Status = case_when(
      has_TRA & has_TRB ~ "TRα & TRβ",
      has_TRA         ~ "TRα Only",
      has_TRB         ~ "TRβ Only",
      TRUE            ~ "None"
    )
  ) %>%
  
  arrange(desc(TCR_Status == "None"))

custom_colors <- c(
  "TRα & TRβ"  = "maroon",
  "TRα Only"   = "orange",
  "TRβ Only"   = "skyblue2",
  "None"       = "grey90"
)

final_plot <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    aes(color = TCR_Status),
    size = 0.025
  ) +
  scale_color_manual(
    values = custom_colors,
    name = "TCR",
    # Only include the desired categories in the legend
    breaks = c("TRα & TRβ", "TRα Only", "TRβ Only")
  ) +
  theme_void() +
  monocle3:::monocle_theme_opts() +
  guides(color = guide_legend(override.aes = list(size = 3)))

print(final_plot)
ggsave("Plots/UMAP_by_TCR.png",dpi=900,height=3,width=3)

cds_hash<-estimate_size_factors(cds_hash)
plot_cells(cds = cds_hash,genes = "IL2",scale_to_range = FALSE,cell_size=.5,label_cell_groups = FALSE) + 
  theme_void() +
  viridis::scale_color_viridis(option = "magma")
plot_cells(cds = cds_hash,genes = "PDCD1",scale_to_range = FALSE,cell_size=.5,label_cell_groups = FALSE) + 
  theme_void() +
  viridis::scale_color_viridis(option = "magma")
plot_cells(cds = cds_hash,genes = "IL2RA",scale_to_range = FALSE,cell_size=.5,label_cell_groups = FALSE) + 
  theme_void() +
  viridis::scale_color_viridis(option = "magma")
plot_cells(cds_hash)+ facet_wrap(~treatment)+geom_density_2d(color = "black") 

plot_cells(cds = cds_hash,genes = c("IL2", "IL2RA", "PDCD1"),scale_to_range = FALSE,cell_size=.5,label_cell_groups = FALSE) + 
  theme_void() +
  viridis::scale_color_viridis(option = "magma")
