plot_percent_cells_positive_pvals <- function (cds_subset, group_cells_by = NULL, min_expr = 0, nrow = NULL, 
                                               ncol = 1, panel_order = NULL, plot_as_count = FALSE, label_by_short_name = TRUE, 
                                               normalize = TRUE, plot_limits = NULL, bootstrap_samples = 100, 
                                               conf_int_alpha = 0.95) 
{
  # 1. Standard Monocle3 setup and assertions for plot_percent_cells_positive
  target_fraction_low <- target_fraction_high <- feature_label <- p_label <- NULL
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  
  marker_exprs <- SingleCellExperiment::counts(cds_subset)
  
  # 2. Create 'marker_exprs_melted'
  if (normalize) {
    marker_exprs <- Matrix::t(Matrix::t(marker_exprs)/size_factors(cds_subset))
    marker_exprs_melted <- reshape2::melt(round(as.matrix(marker_exprs), digits = 4))
  } else {
    marker_exprs_melted <- reshape2::melt(as.matrix(marker_exprs))
  }
  
  colnames(marker_exprs_melted) <- c("f_id", "cell", "expression")
  c_data <- as.data.frame(colData(cds_subset))
  r_data <- as.data.frame(rowData(cds_subset))
  c_data$cell_id_internal <- rownames(c_data)
  r_data$feat_id_internal <- rownames(r_data)
  
  marker_exprs_melted <- merge(marker_exprs_melted, c_data, by.x = "cell", by.y = "cell_id_internal")
  marker_exprs_melted <- merge(marker_exprs_melted, r_data, by.x = "f_id", by.y = "feat_id_internal")
  
  # 3. Gene labeling
  if (label_by_short_name && !is.null(marker_exprs_melted$gene_short_name)) {
    marker_exprs_melted$feature_label <- marker_exprs_melted$gene_short_name
    marker_exprs_melted$feature_label[is.na(marker_exprs_melted$feature_label)] <- marker_exprs_melted$f_id
  } else {
    marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
  }
  
  if (!is.null(panel_order)) {
    marker_exprs_melted$feature_label <- factor(marker_exprs_melted$feature_label, levels = panel_order)
  }
  
  if (is.null(group_cells_by)) {
    stop("group_cells_by must be provided for a between-group comparison.")
  }
  
  # 4. Between-group Bootstrap Logic
  # Split by gene only to keep groups together for comparison
  gene_list <- split(marker_exprs_melted, marker_exprs_melted$feature_label, drop = TRUE)
  
  marker_counts <- lapply(gene_list, function(sub_data) {
    if(nrow(sub_data) == 0) return(NULL)
    
    lvls <- unique(as.character(sub_data[[group_cells_by]]))
    if(length(lvls) != 2) {
      stop(paste("Gene", sub_data$feature_label[1], "does not have exactly 2 groups."))
    }
    
    # Bootstrap the DIFFERENCE between groups
    # Statistic: Fraction(Group 2) - Fraction(Group 1)
    b_obj_diff <- boot::boot(data = sub_data, statistic = function(d, i) {
      d_sub <- d[i, ]
      m1 <- mean(d_sub$expression[d_sub[[group_cells_by]] == lvls[1]] > min_expr)
      m2 <- mean(d_sub$expression[d_sub[[group_cells_by]] == lvls[2]] > min_expr)
      return(m2 - m1)
    }, R = bootstrap_samples)
    
    # Calculate P-value for the difference (Null: difference = 0)
    p_val <- boot.pval::boot.pval(b_obj_diff, type = "perc", theta_null = 0)
    
    # Calculate individual group stats for bar heights and error bars
    group_stats <- do.call(rbind, lapply(lvls, function(l) {
      grp_data <- sub_data[sub_data[[group_cells_by]] == l, ]
      b_grp <- boot::boot(data = grp_data, statistic = function(d, i) {
        mean(d[i, ]$expression > min_expr)
      }, R = bootstrap_samples)
      
      data.frame(
        feature_label = sub_data$feature_label[1],
        group_var = l,
        target_fraction_mean = b_grp$t0,
        target_fraction_low = stats::quantile(b_grp$t, (1 - conf_int_alpha)/2),
        target_fraction_high = stats::quantile(b_grp$t, 1 - (1 - conf_int_alpha)/2),
        p_value = p_val,
        stringsAsFactors = FALSE
      )
    }))
    return(group_stats)
  }) %>% do.call(rbind, .)
  
  colnames(marker_counts)[which(names(marker_counts) == "group_var")] <- group_cells_by
  
  # 5. Format P-value labels (Small text, 2 sig figs)
  marker_counts$p_label <- sapply(marker_counts$p_value, function(p) {
    p_val_fmt <- format(signif(p, 2))
    if (p < 0.05) return(paste0(p_val_fmt, "*"))
    return(p_val_fmt)
  })
  
  # 6. Final Plot construction
  marker_counts$target_fraction_mean <- marker_counts$target_fraction_mean * 100
  marker_counts$target_fraction_low <- marker_counts$target_fraction_low * 100
  marker_counts$target_fraction_high <- marker_counts$target_fraction_high * 100
  
  # Ensure the group variable is a factor for proper x-axis ordering
  marker_counts[[group_cells_by]] <- factor(marker_counts[[group_cells_by]])
  
  qp <- ggplot(marker_counts, aes_string(x = group_cells_by, y = "target_fraction_mean", fill = group_cells_by)) +
    facet_wrap(~feature_label, nrow = nrow, ncol = ncol, scales = "free_y") +
    geom_bar(stat = "identity") + 
    geom_linerange(aes(ymin = target_fraction_low, ymax = target_fraction_high)) + 
    geom_text(aes(label = ifelse(duplicated(feature_label), p_label, ""), y = Inf), 
              vjust = 1.5, size = 1.5) +
    ylab("Cells (percent)") +
    monocle3:::monocle_theme_opts()
  
  return(qp)
}
