# Visualize sequencing quality metrics (duplication rate and read distribution).

library(dplyr)
library(ggplot2)
library(monocle3)

# Define input and output directories
input_dir <- "path/to/input" # Change this to your input directory
output_dir <- "path/to/output" # Change this to your output directory
# setwd(output_dir) # Optional

# --------------------------------------------------------------------------------
# Duplication Rate vs Reads per Cell
# --------------------------------------------------------------------------------
df_all <- readRDS(file.path(input_dir, "df_all.RDS"))

df_all$method <- sapply(df_all$dataset, function(x){
  ifelse(x == "DT", "dT", ifelse(x == "FF", "TRTL", "TRTL-ON"))
})

ggplot(df_all, aes(x = reads_per_cell, y = dup_rate, color = method)) +
  geom_point(size = 1.5) +
  geom_line(size = .5, aes(group = method)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    values = c(
      "dT" = "dimgrey",   
      "TRTL" = "navy",   
      "TRTL-ON" = "deepskyblue2"  
    )
  ) +
  labs(
    x = "Mean reads per cell",
    y = "Duplication rate",
    color = "Method"
  ) +
  monocle3:::monocle_theme_opts() +
  theme(text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.5,"line")) +
  guides(colour = guide_legend(override.aes = list(size=1)))

ggsave(file.path(output_dir, "dup_rate_vs_reads_per_cell.png"), width = 2.5,height = 1.5,dpi = 900)

# --------------------------------------------------------------------------------
# Read Distribution (Bar Plot)
# --------------------------------------------------------------------------------
df_long <- readRDS(file.path(input_dir, "df_long.RDS"))

df_long$method <- sapply(df_long$method, function(x){ifelse(x == "targeted", "TRTL", x)})
colnames(df_long) <- c("method", "Category", "total", "proportion")

df_long$Category <- sapply(df_long$Category, function(x){
  if(x == "MultiMappingReads") return("Multi Mapped")
  if(x == "UniqueMappedReads") return("Uniquely Mapped")
  if(x == "UnmappedReadsTooShort") return("Unmapped")
  return(NA)
  })
                                                              
color_scheme <- c(
  "Uniquely Mapped" = "dimgrey",  
  "Multi Mapped" = "grey90",  
  "Unmapped" = "grey20"  
)

ggplot(df_long, aes(x = method, y = proportion, fill = Category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +  
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Method", y = "Fraction reads\n(%)") +
  theme_minimal() +
  scale_fill_manual(values = color_scheme) +  
  theme(
    legend.position = "right", 
    legend.direction = "vertical", 
    text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.width = unit(0.5,"line"),
    legend.key.height = unit(0.5,"line"),
    panel.grid = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size=1)))

ggsave(file.path(output_dir, "Read_distribution.png"), dpi = 900, height = 1.5, width = 2.5)
