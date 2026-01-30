group_cell_subtypes <- function(colData_df) {
  #' Aggregates cell subtypes based on the last word of the class_name.
  #'
  #' This function takes a cell metadata data frame, extracts the core cell
  #' type from the 'class_name' column (the word after the last space),
  #' and then groups these types into major categories (Glutamatergic,
  #' GABAergic, Other).
  #'
  #' @param colData_df A data frame containing cell metadata (e.g.,
  #'   as.data.frame(colData(cds))). Must contain a column named 'class_name'.
  #' @return A data frame with two new columns: 'simplified_cell_type' (the
  #'   last word of class_name) and 'major_group_cell_type' (Glutamatergic,
  #'   GABAergic, or Other).
  
  # 1. Extract the last word from the class_name (everything after the last space)
  # The pattern ".*\\s" matches any characters followed by the last space.
  # str_replace replaces the matched pattern with an empty string, leaving only
  # the last word (the subtype).
  
  df_grouped <- colData_df %>%
    mutate(
      simplified_cell_type = str_replace(class_name, ".*\\s", "")
    )
  
  df_grouped <- df_grouped %>%
    mutate(
      major_group_cell_type = case_when(
        str_detect(simplified_cell_type, "Glut") ~ "Glutamatergic",
        str_detect(simplified_cell_type, "GABA") ~ "GABAergic",
        TRUE ~ simplified_cell_type 
      )
    ) %>%
    
    select(
      starts_with("class_"), 
      class_name,
      simplified_cell_type, 
      major_group_cell_type
    )
  
  return(df_grouped)
}
