library(tidyverse)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)



args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Please provide input_dir and output_dir as arguments.")
}
input_dir <- args[1]
output_dir <- args[2]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

files <- list.files(path = input_dir, pattern = "*.tsv", full.names = TRUE)

# Read all files and add the PCR column
data_list <- files %>%
  map(~ read_tsv(.)) %>%
  map2(files, ~ mutate(.x, PCR = tools::file_path_sans_ext(basename(.y))))
filtered_data_list <- discard(data_list, ~ nrow(.) == 0)
# Combine all data frames into one
combined_data <- bind_rows(filtered_data_list)

# View the combined data
head(combined_data)
combined_data$cell <-paste(combined_data$tagValueCELL2ID, combined_data$tagValueCELL1ID, substr(combined_data$PCR, 12, 14), sep = "_")
length(unique(combined_data$cell)) 
unique(colnames(combined_data))

saveRDS(combined_data, file.path(output_dir, "TCR_wide_format_mixcr.RDS"))
filtered_data <- combined_data %>%
  group_by(PCR) %>%
  mutate(pcr_total_reads = n()) %>%
  group_by(PCR, targetSequences) %>%
  mutate(
    target_seq_reads = n(),
    percentage = (target_seq_reads / pcr_total_reads) * 100
  ) %>%
  ungroup() %>%
  filter(percentage >= 2) %>%
  select(-pcr_total_reads, -target_seq_reads, -percentage)

TCR_summary <- filtered_data %>% #combined_data %>%
  group_by(cell, topChains) %>%  # Group by cell and chain type (TRA/TRB)
  slice_max(order_by = readFraction, n = 1, with_ties = FALSE) %>%  # Keep highest readFraction row per chain
  select(cell, topChains, targetSequences, uniqueMoleculeCount) %>%  # Keep relevant columns
  pivot_wider(
    names_from = topChains, 
    values_from = c(targetSequences, uniqueMoleculeCount)
  ) %>%
  ungroup()

TCR_summary <- TCR_summary %>%
  rename_with(~ gsub("targetSequences_", "", .), starts_with("targetSequences")) %>%
  rename_with(~ gsub("uniqueMoleculeCount_", "UMI_", .), starts_with("uniqueMoleculeCount")) 

saveRDS(TCR_summary, file.path(output_dir, "TCR_summary.RDS"))

