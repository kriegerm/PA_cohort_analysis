# ---- Setup ----
set.seed(1234)
date = Sys.Date()
 
#Make output directory for data saved from functions
# Specify the directory path
dir_path <- "../saved_analysis_files"

# Check if the directory exists
if (!dir.exists(dir_path)) {
  # Create the directory if it doesn't exist
  dir.create(dir_path, recursive = TRUE)  # 'recursive' allows creating parent directories if they don't exist
  message("Directory created: ", dir_path)
} else {
  message("Directory already exists: ", dir_path)
}



# ---- Preprocessing Functions ----

## Construct phyloseq function

#Example:
#Define variables and file locations
#biom_location <- "PRJNA822685_OSCC/Qiime2/blast/Files_For_Phyloseq/feature_table_w_taxonomy.biom"
#tree_location <- "PRJNA822685_OSCC/Qiime2/Tree/Unfiltered_Rooted_tree_for_phyloseq/tree.nwk"
#sampledata_location <- "PRJNA822685_OSCC/Qiime2/Files_For_Phyloseq/metadata_phyloseq.tsv"
#construct_phyloseq("phy_obj_oscc", biom_location, tree_location, sampledata_location)

construct_phyloseq <- function(phyloseq_object_name, biom_location, tree_location, sampledata_location) {
  # Validate inputs
  if (!file.exists(biom_location)) stop("The biom file does not exist.")
  if (!file.exists(tree_location)) stop("The tree file does not exist.")
  if (!file.exists(sampledata_location)) stop("The sample data file does not exist.")
  
  ps1 <- import_biom(biom_location) 
  tree <- read_tree(tree_location)
  
  #Add the Taxonomic Ranks to the different Taxonomic levels
  colnames(tax_table(ps1)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  #Combine all F. nucleatum subspecies
  tax_tab_extracted <- tax_table(ps1) %>% as.matrix() %>% as.data.frame() %>%
    mutate(Species = str_replace(Species, "^s__nucleatum_subsp\\._.*", "s__nucleatum"))
  tax_tab <- tax_table(as.matrix(tax_tab_extracted))
  
  OTU_tab <- ps1@otu_table 
  
  #Generate the Phyloseq Object
  phylo_obj <- phyloseq(OTU_tab, tax_tab)
  sample_data(phylo_obj) <- import_qiime_sample_data(sampledata_location)
  
  #Update the tax table so that genus and species are joined for the species column
  tax_df <- phylo_obj@tax_table %>% as.matrix() %>% as.data.frame() %>%
        unite(Genus_species, c(Genus, Species), sep = "_", remove = FALSE) %>%
        dplyr::select(-c("Species")) %>%
        dplyr::rename("Species" = "Genus_species") %>%
        dplyr::select(-"Species", "Species")
  
  # Ensure row names are preserved and match the original tax table
  rownames(tax_df) <- taxa_names(phylo_obj)
  # Convert back to matrix
  new_tax_matrix <- as.matrix(tax_df)
  # Create a new tax_table object
  new_tax_table <- tax_table(new_tax_matrix)
  # Ensure that taxa names are consistent
  taxa_names(new_tax_table) <- taxa_names(phylo_obj)
  # Update the phyloseq object with the new tax table
  phylo_obj@tax_table <- new_tax_table
  
  # Assign the object to the global environment to output
  assign(phyloseq_object_name, phylo_obj, envir = .GlobalEnv)
}

## ---- Filter Unmatched----
#The metadata in the Pediatric and Adult Abscess study differ slightly, so I just made two different functions to retain ONLY samples with 1 plaque and 1 abscess each.

## Pediatric Study: Filter for unmatched samples
# This function is used inside the function filter_phyloseq
filter_unmatched_samples_Pediatric <- function(phyloseq_obj) {
  # Extract the sample data
  sample_data_df <- sample_data(phyloseq_obj) %>% as.matrix() %>%
    as.data.frame()
  
  # Ensure the necessary columns exist
  required_columns <- c("Sample", "Type")
  if (!all(required_columns %in% colnames(sample_data_df))) {
    stop("The sample data must contain 'Sample' and 'Type' columns.")
  }
  
  # Extract numeric parts of sample names
  numeric_sample_names <- gsub("-.*", "", sample_names(phyloseq_obj))
  sample_data_df$NumericSample <- gsub("-.*", "", sample_data_df$Sample)

  # Count the number of unique 'Type' entries for each numeric sample
  type_counts <- sample_data_df %>%
    dplyr::select(c("NumericSample", "Type")) %>%
    dplyr::group_by(NumericSample) %>%
    dplyr::summarize(num_types = n_distinct(Type), .groups = "drop")
  
  # Filter for numeric samples with only one type
  single_type_samples <- type_counts %>%
    filter(num_types == 1) %>%
    pull(NumericSample)
  
  print(single_type_samples)
  
  # Filter phyloseq object by numeric samples
  filtered_phyloseq_obj <- prune_samples(
    !(numeric_sample_names %in% single_type_samples),
    phyloseq_obj
  )
  
  return(filtered_phyloseq_obj)
}


## ---- Filter Phyloseq ----
#Example:
#filter_phyloseq(phy_obj_oscc, "Pediatric", Contam_g, Contam_f, Contam_s)

filter_phyloseq <- function(phylo_obj, study, Contam_g, Contam_f, Contam_s) {
  # Capture the name of the input variable
  phylo_name <- deparse(substitute(phylo_obj))
  
    # Ensure the 'Study' column exists in the sample data
  if (!"Study" %in% colnames(sample_data(phylo_obj))) {
    stop("The 'Study' column is not present in the sample data of the phyloseq object.")
  }
  
  #Filter out samples only in the study - important if you have other studies included in the metadata
  phylo_obj_f <- switch(study,
                      "OSCC" = subset_samples(phylo_obj, Study == "OSCC"),
                      "Pediatric" = subset_samples(phylo_obj, Study == "Pediatric"),
                      "AdultAbscess" = subset_samples(phylo_obj, Study == "AdultAbscess"))
    

  #Filter Unassigned
  phylo_obj_f <- subset_taxa(phylo_obj_f, Phylum != "Unassigned")

  #filter zeros
  phylo_obj_f <- prune_taxa(taxa_sums(phylo_obj_f) > 0, phylo_obj_f)
  
  #fix taxa with Microviz
  phylo_obj_f <- tax_fix(phylo_obj_f)
  
  #Remove values with less than ___ total reads:
  phylo_obj_f<- prune_samples(sample_sums(phylo_obj_f)>= 1000, phylo_obj_f)
  
   #Remove all samples that are not paired (i.e. only have one plaque and one abscess)
  if (study == "Pediatric") {
    phylo_obj_f <- filter_unmatched_samples_Pediatric(phylo_obj_f)
  }else if (study != "Pediatric") {
    phylo_obj_f <- filter_unmatched_samples(phylo_obj_f)
  }

  #Remove contaminates
  phylo_obj_f <- subset_taxa(phylo_obj_f, !(Genus %in% Contam_g)) 
  phylo_obj_f <- subset_taxa(phylo_obj_f, !(Family %in% Contam_f)) 
  phylo_obj_f <- subset_taxa(phylo_obj_f, !(Species %in% Contam_s))
    
  #Output
  assign(paste0(phylo_name, "_f"), phylo_obj_f, envir = .GlobalEnv)
}


## ---- Process Phyloseq ----


process_phyloseq <- function(phylo_obj_fs) {
  
  # Capture the name of the input variable
  phylo_name <- deparse(substitute(phylo_obj_fs))
  
  # Transform to relative abundance with translated pseudo count
  phylo_obj_fs_t <- transform_sample_counts(phylo_obj_fs, function(x) (x + 0.000001 - min(x)))
  phylo_obj_fs_t_n <- transform_sample_counts(phylo_obj_fs_t, function(x) (x) / sum(x))
  
  # Taxonomic agglomeration at different levels
  phylo_obj_fs_t_n_pglom <- tax_glom(phylo_obj_fs_t_n, taxrank = "Phylum", NArm = TRUE) 
  phylo_obj_fs_t_n_gglom <- tax_glom(phylo_obj_fs_t_n, taxrank = "Genus", NArm = TRUE)
  phylo_obj_fs_t_n_sglom <- tax_glom(phylo_obj_fs_t_n, taxrank = "Species", NArm = TRUE)

  # Convert to data.frame
  phylo_obj_fs_t_n_pglom_df <- psmelt(phylo_obj_fs_t_n_pglom) 
  phylo_obj_fs_t_n_gglom_df <- psmelt(phylo_obj_fs_t_n_gglom)
  phylo_obj_fs_t_n_sglom_df <- psmelt(phylo_obj_fs_t_n_sglom)
  
  # Prepare for CSV export
  species_csv <- phylo_obj_fs_t_n_sglom_df %>%
    dplyr::select(c("Sample", "Abundance", "Species", "Type")) %>%
    pivot_wider(names_from = "Species", values_from = "Abundance")
  
  genus_csv <- phylo_obj_fs_t_n_gglom_df %>%
    dplyr::select(c("Sample", "Abundance", "Genus", "Type")) %>%
    pivot_wider(names_from = "Genus", values_from = "Abundance")
  
  phylum_csv <- phylo_obj_fs_t_n_pglom_df %>%
    dplyr::select(c("Sample", "Abundance", "Phylum", "Type")) %>%
    pivot_wider(names_from = "Phylum", values_from = "Abundance")
  
  
  # Dynamically name and assign outputs
  assign(paste0(phylo_name, "_t_n"), phylo_obj_fs_t_n, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_pglom_df"), phylo_obj_fs_t_n_pglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_gglom_df"), phylo_obj_fs_t_n_gglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_sglom_df"), phylo_obj_fs_t_n_sglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_species_csv"), species_csv, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_genus_csv"), genus_csv, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_phylum_csv"), phylum_csv, envir = .GlobalEnv)
  
}

# ---- Basic Taxa Plots ----

#plot_all_taxa(ps_obj_ped_fs, "TSS", "condition", "Genus", "dot", plot_colors)
plot_all_taxa <- function(ps_obj,
                          transformation = "TSS",
                          group,
                          taxa_level,
                          plot_type = "box",
                          plot_colors = NULL,
                          taxa_filter = NULL) {
  
  if (transformation == "counts") {
  plot_data <- tax_glom(ps_obj, taxa_level) %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
    mutate(Count_Type = "Count") 

  }else if (transformation == "TSS") {
  plot_data <- tax_glom(ps_obj, taxa_level) %>%
    transform_sample_counts(., function(x) (x + 0.000001 - min(x))) %>%
    transform_sample_counts(., function(x) (x) / sum(x)) %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
    mutate(Count_Type = "RelAb")
  
  }else if (transformation == "CLR") {
  plot_data<- tax_glom(ps_obj, taxa_level) %>%
    microbiome::transform(., transform = "clr", target="sample") %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
    mutate(Count_Type = "CLR")    

  } else {
      stop("Invalid transformation type. Please specify 'TSS', 'CLR', or 'counts'.")
  }
  
  # Filter for selected taxa, if provided
  if (!is.null(taxa_filter)) {
    plot_data <- plot_data %>%
      filter(!!sym(taxa_level) %in% taxa_filter) %>%
      mutate(!!sym(taxa_level) := factor(!!sym(taxa_level), levels = taxa_filter))
  }
  
   # Create dot plot
  if (plot_type == "dot") {
  plot <- plot_data %>%
    ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    facet_wrap(vars(!!sym(taxa_level), Count_Type), ncol=6, scales = "free") + # Use facet_wrap with free scales
    theme_bw(base_size = 10) +
    scale_color_manual(values = plot_colors) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      
    )
  
   #Box plot
  }else if (plot_type == "box") {
    plot <- plot_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_boxplot(alpha = 0.5, outlier.shape = 8) + # Add boxplot
      #geom_jitter(width = 0.2, alpha = 0.7, size = 1.5) + # Add jittered points
      facet_wrap(vars(!!sym(taxa_level), Count_Type), ncol = 6, scales = "free") + # Facet
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
  #Violin
  }else if (plot_type == "violin") {
    plot <- plot_data %>%
        ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
        geom_violin(alpha = 0.5, scale = "width") + # Add violin plot with partial transparency
        geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 1.5, alpha = 0.7) + # Add jittered points
        facet_wrap(vars(!!sym(taxa_level), Count_Type), ncol = 6, scales = "free") + # Facet
        theme_bw(base_size = 10) +
        scale_color_manual(values = plot_colors) +
        scale_fill_manual(values = plot_colors) +
        theme(
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
  } else {
      stop("Invalid plot type. Please specify 'box', 'violin', 'density' or 'dot'.")
  }
  return(list(plot = plot))
}



# ---- Taxa Trends ----

# Plot the taxa trends with bubbles
plot_taxa_shifts <- function(top_number = 10,
                             input_df = ps_fs_genus_csv,
                             bubble_max_size = 10,
                             count_max = 25) {
  
  
  # --------- Long form + per-taxon stats
  df_long <- input_df %>%
    dplyr::select(-Sample) %>%
    tidyr::pivot_longer(-Type, names_to = "Taxa", values_to = "Abundance")
  
  stopifnot(all(c("Type","Taxa","Abundance") %in% names(df_long)))
  
  summ <- df_long %>%
    dplyr::group_by(Taxa) %>%
    dplyr::summarise(
      n_plaque      = sum(Type == "Plaque"),
      n_abscess     = sum(Type == "Abscess"),
      median_plaque = median(Abundance[Type == "Plaque"],  na.rm = TRUE),
      median_abscess= median(Abundance[Type == "Abscess"], na.rm = TRUE),
      mean_plaque   = mean(Abundance[Type == "Plaque"],    na.rm = TRUE),
      mean_abscess  = mean(Abundance[Type == "Abscess"],   na.rm = TRUE),
      p = tryCatch(
        wilcox.test(Abundance[Type == "Plaque"], Abundance[Type == "Abscess"],
                    exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      diff_median = median_plaque - median_abscess,
      diff_mean   = mean_plaque   - mean_abscess,
      direction   = dplyr::case_when(
        diff_median >  0 ~ "Plaque",
        diff_median <  0 ~ "Abscess",
        TRUE             ~ "Tie"
      ),
      padj = p.adjust(p, method = "BH")
    )
  
  # ----- Pick top N per direction + clean taxon labels
  p_to_stars <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", ""))))
  }
  
  topN <- summ %>%
    dplyr::filter(padj <= .05, direction %in% c("Plaque","Abscess")) %>%
    dplyr::group_by(direction) %>%
    dplyr::slice_max(order_by = abs(diff_median), n = top_number, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Taxa_clean = Taxa %>%
        stringr::str_remove("^g__") %>%
        stringr::str_remove_all("s__") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_replace_all("G ", "G") %>%
        stringr::str_remove_all("\\[|\\]"),
      signed_effect = diff_median,
      star = p_to_stars(padj)
    )
  
  # y order (smallest signed_effect at bottom, largest at top after coord_flip)
  y_levels <- topN %>%
    dplyr::arrange(signed_effect) %>%
    dplyr::pull(Taxa_clean)
  y_levels_plot <- rev(y_levels)  # this is the visual top→bottom
  
  # offset for star positioning
  offset <- 0.1 * max(abs(topN$signed_effect), na.rm = TRUE)
  topN <- topN %>%
    dplyr::mutate(
      Taxa_clean = factor(Taxa_clean, levels = y_levels),   # base order
      label_y = signed_effect + ifelse(signed_effect >= 0, offset, -offset)
    )
  
  # ------ Bar plot (left). Uses y_levels_plot for axis.
  # define palette for bars
  colors_all <- c(Plaque = "#3185FC", Abscess = "#FF495C", Tie = "grey70")
  
  bar_plot <- ggplot(
    topN,
    aes(x = Taxa_clean, y = signed_effect, fill = direction)
  ) +
    geom_col(color = "black") +
    geom_text(aes(y = label_y, label = star), size = 5, vjust = .8) +
    coord_flip() +
    # vertical gridlines across categories (optional aesthetic)
    geom_vline(xintercept = seq(1.5, length(y_levels) - 0.5, 1),
               color = "black", size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(limits = y_levels_plot) +  # enforce exact shown order
    labs(x = NULL, y = "", fill = "Enriched in") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "italic", color="black", size=14),
      axis.text.x = element_text(angle = 45, color="black", size=11, hjust=1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = colors_all)
  
  # ---------- Build bubble counts (A vs P) with SAME y
  tax_pull <- topN %>% dplyr::pull(Taxa)
  
  data_pull <- input_df %>%
    tidyr::pivot_longer(cols = starts_with("g_"), names_to = "Taxa", values_to = "relab") %>%
    dplyr::filter(Taxa %in% tax_pull) %>%
    dplyr::mutate(
      Taxa_clean = Taxa %>%
        stringr::str_remove("^g__") %>%
        stringr::str_remove_all("s__") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_remove_all("\\[|\\]"),
      Sample = stringr::str_replace_all(Sample, "-P|-[Aa]", "")
    ) %>%
    tidyr::pivot_wider(names_from = "Type", values_from = "relab") %>%
    dplyr::mutate(
      Ratio     = Abscess / Plaque,
      RatioFlag = dplyr::if_else(Ratio > 1, 1, -1, missing = NA_integer_)
    )
  
  df_summary <- data_pull %>%
    dplyr::group_by(Taxa_clean) %>%
    dplyr::summarise(
      Abscess_enriched = sum(RatioFlag == 1,  na.rm = TRUE),
      Plaque_enriched  = sum(RatioFlag == -1, na.rm = TRUE),
      .groups = "drop"
    )
  
  # keep topN order only (and ensure exact order)
  df_summary <- df_summary %>%
    dplyr::slice(match(y_levels, Taxa_clean))
  
  # long for plotting
  bubble_df <- df_summary %>%
    tidyr::pivot_longer(
      cols = c(Abscess_enriched, Plaque_enriched),
      names_to = "group2", values_to = "count"
    ) %>%
    dplyr::mutate(
      group2     = dplyr::recode(group2, "Abscess_enriched" = "A", "Plaque_enriched" = "P"),
      Taxa_clean = factor(Taxa_clean, levels = y_levels) # same base order
    )
  
  # scaffold to guarantee both A and P exist for every taxon (keeps row heights identical)
  scaffold <- expand.grid(
    Taxa_clean = factor(y_levels, levels = y_levels),
    group2     = c("A","P"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  
  bubble_df_full <- scaffold %>%
    dplyr::left_join(bubble_df, by = c("Taxa_clean","group2")) %>%
    dplyr::mutate(count = tidyr::replace_na(count, 0L))
  
  bubble_A <- bubble_df_full %>% dplyr::filter(group2 == "A")
  bubble_P <- bubble_df_full %>% dplyr::filter(group2 == "P")
  
  # -------- Bubble plot (right) with the SAME y order
  bubble_plot <- ggplot() +
    # Abscess (red)
    geom_point(
      data = bubble_A,
      aes(x = group2, y = Taxa_clean, size = count, fill = count),
      shape = 21, alpha = 0.9, stroke = 0.4, color = "black"
    ) +
    scale_fill_gradientn(
      colors = c("white", "#FF495C"),
      limits = c(0, count_max),
      oob = scales::squish,
      guide = "none" ,                 # set to colorbar to show a legend
      name = "Abscess (A)"
    ) +
    ggnewscale::new_scale("fill") +
    # Plaque (blue)
    geom_point(
      data = bubble_P,
      aes(x = group2, y = Taxa_clean, size = count, fill = count),
      shape = 21, alpha = 0.9, stroke = 0.4, color = "black"
    ) +
    scale_fill_gradientn(
      colors = c("white", "#3185FC"),
      limits = c(0, count_max),
      oob = scales::squish,
      guide = "none",                  # set to colorbar to show a legend
      name = "Plaque (P)"
    ) +
    scale_size_area(
      name = "Count",
      max_size = bubble_max_size,
      limits = c(0, count_max)
    ) +
    scale_x_discrete() +
    # give a hair of vertical breathing room so edge bubbles don't touch frame
    scale_y_discrete(limits = y_levels_plot, drop = TRUE, expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.border = element_blank(),   # <-- remove the panel border
      axis.line    = element_blank(),   # (optional) remove axis lines too
      text = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(face = "italic", size = 12),
      legend.position = "right",
      legend.box = "horizontal"
    )
  
  # ------- Stitch (bar left, bubbles right)
  combined <- bar_plot +
    (bubble_plot +
       labs(x = NULL) +
       theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank())) +
    plot_layout(widths = c(3, 1)) 
  
  list(
    bar_plot     = bar_plot,
    bubble_plot  = bubble_plot,
    combined_plot= combined,
    data         = summ,
    bubble_data  = bubble_df_full
  )
}


# ---- Ordinations ----

## ---- Compare methods ----

#BIG function to compare different ordination methods 
compare_ordination_methods <- function(phyloseq_obj,
                                       rank_transformations = c("Species", "Genus", "Phylum"),
                                       trans_types = c("relab", "clr", "identity", "log"),
                                       dist_cal_types = c("euclidean", "bray", "jaccard"),
                                       ord_calc_methods = c("PCA", "PCoA", "NMDS"),
                                       No_axes_to_correlated_with_LibSize = 2,
                                       variance_threshold = 0.80) { #find the number of axes needed for PCA and PCoA that explain this percent of variance  
  
  #Create a parameter hash so that you can just load results if they've already been run in an identical way
  param_hash <- digest(list(
    phyloseq_obj = phyloseq_obj,
    rank_transformations = rank_transformations,
    trans_types = trans_types,
    dist_cal_types = dist_cal_types,
    ord_calc_methods = ord_calc_methods,
    No_axes_to_correlated_with_LibSize = No_axes_to_correlated_with_LibSize,
    variance_threshold = variance_threshold  # Include in hash for reproducibility
  ))
  
  result_file <- file.path("../saved_analysis_files/", paste0("ordination_parameter_sweep_result_", param_hash, ".rds"))
  
  if (file.exists(result_file)) {
    message("Analysis already run. Loading results...")
    results_df <- readRDS(result_file)
    return(list(results_df = results_df))
    
  } else {
    message("Running analysis...")
  }
  
  #Define the CLR transformation
  clr_transform <- function(mat, pseudocount = 1) {
    mat <- mat + pseudocount
    log_mat <- log(mat)
    gm <- rowMeans(log_mat)
    sweep(log_mat, 1, gm)
  }
  
  validate_transformation_distance <- function(trans_type, dist_type, ord_method) {
    ord_method <- toupper(ord_method)
    trans_type <- tolower(trans_type)
    dist_type <- tolower(dist_type)
    
    ## Warnings and strict incompatibilities
    
    # CLR requires Euclidean
    if (trans_type == "clr" && dist_type != "euclidean") {
      stop("CLR transformation is only compatible with Euclidean distance.")
    }
    
    # NMDS is discouraged with Euclidean 
    if (ord_method == "NMDS" && dist_type == "euclidean") {
      warning("NMDS with Euclidean distance is discouraged in ecological/microbiome data. Prefer Bray or Jaccard.")
    }
    
    # NMDS with CLR is not supported (theoretically questionable)
    if (ord_method == "NMDS" && trans_type == "clr") {
      stop("NMDS is not compatible with CLR-transformed data.")
    }
    
    # Identity (raw counts) with Bray is discouraged
    if (trans_type == "identity" && dist_type == "bray") {
      stop("Using raw counts with Bray-Curtis distance may give misleading results. Consider using relative abundance or log transformation.")
    }

    # Log transformation with Jaccard is questionable
    if (trans_type == "log" && dist_type == "jaccard") {
      stop("Jaccard distance with log-transformed data is rarely appropriate. Use presence/absence or relative abundance instead.")
    }
    
    # Log transformation with Bray: acceptable, but might want a warning
    if (trans_type == "log" && dist_type == "bray") {
      stop("Log-transformed data with Bray-Curtis distance is valid, but ensure zero handling is appropriate.")
    }
    
    return(dist_type)
  }
  
  #Initiate empty results list and set counter
  results_list <- list()
  counter <- 1
  
  #Begin looping
  for (rank_trans in rank_transformations) {
    for (trans_type in trans_types) {
      for (ord_method in ord_calc_methods) {
        ord_method_upper <- toupper(ord_method)
        
        #PCA only uses euclidean, so if  necessary set that
        dist_types_to_use <- if (ord_method_upper == "PCA") "euclidean" else dist_cal_types
        
        for (dist_type in dist_types_to_use) {
          message(sprintf("Trying: %s | %s | %s | %s", rank_trans, trans_type, dist_type, ord_method))
          
          #Try the combination - if it's not compatible that's ok, it'll skip it
          try({
            phylo_trans <- tax_glom(phyloseq_obj, taxrank = rank_trans)
            otu_mat <- as.matrix(otu_table(phylo_trans))
            if (taxa_are_rows(phylo_trans)) {
              otu_mat <- t(otu_mat)
            }
            
            if (trans_type == "identity") {
              otu_trans <- otu_mat
            } else if (trans_type == "log") {
              otu_trans <- log1p(otu_mat)
            } else if (trans_type == "relab") {
              otu_trans <- sweep(otu_mat, 1, rowSums(otu_mat), FUN = "/")
            } else if (trans_type == "clr") {
              otu_trans <- clr_transform(otu_mat)
            } else {
              stop("Unsupported transformation.")
            }
            
            #Remove rows where the variance is 0
            otu_trans <- otu_trans[, apply(otu_trans, 2, var) != 0]
            
            #Use the function defined at the beginning to validate that the transformation and distance type will work
            dist_type_valid <- validate_transformation_distance(trans_type, dist_type, ord_method)
            
            #Intiate empty components to fill with the functions below
            scores_df <- NULL
            var_expl <- NULL
            stress <- NA
            n_components_for_variance <- NA 
            
            if (ord_method_upper == "PCA") {
              ord <- prcomp(otu_trans, center = TRUE, scale. = TRUE)
              scores_df <- as.data.frame(ord$x)
              colnames(scores_df) <- paste0("Axis", 1:ncol(scores_df))
              var_expl <- ord$sdev^2 / sum(ord$sdev^2)
              
              # Calculate components needed to explain the variance threshold
              cumulative_var <- cumsum(var_expl)
              n_components_for_variance <- as.numeric(which(cumulative_var >= variance_threshold)[1])
              
            } else if (ord_method_upper == "PCOA") {
              
              #Binarize for Jaccard
              if (dist_type_valid == "jaccard") {
                dist_mat <- vegdist(otu_trans, method = "jaccard", binary=TRUE)
              } else if (dist_type_valid != "jaccard") {
                dist_mat <- vegdist(otu_trans, method = dist_type_valid)
              }
              
              n_samples <- attr(dist_mat, "Size") 
              ord <- cmdscale(dist_mat, k = n_samples - 1, eig = TRUE)
              scores_df <- as.data.frame(ord$points)
              colnames(scores_df) <- paste0("Axis", 1:ncol(scores_df))
              var_expl <- ord$eig / sum(ord$eig)
              
              # Calculate components needed to explain the variance threshold
              cumulative_var <- cumsum(var_expl)
              n_components_for_variance <- as.numeric(which(cumulative_var >= variance_threshold)[1])
              
            } else if (ord_method_upper == "NMDS") {
              
              #Binarize for Jaccard
              if (dist_type_valid == "jaccard") {
                dist_mat <- vegdist(otu_trans, method = "jaccard", binary=TRUE)
              } else if (dist_type_valid != "jaccard") {
                dist_mat <- vegdist(otu_trans, method = dist_type_valid)
              }
              
              ord <- metaMDS(dist_mat, trymax = 100)
              scores_df <- as.data.frame(ord$points)
              colnames(scores_df) <- paste0("Axis", 1:ncol(scores_df))
              stress <- ord$stress
              n_components_for_variance <- as.numeric(2) #Automatically set because there are only two axes - need this to be a value for later
            }
            
            #Join the metadata with the scores_df
            scores_df$SampleID <- rownames(scores_df)
            meta_df <- sample_data(phylo_trans) %>% as.matrix() %>% as.data.frame()
            scores_df <- dplyr::left_join(scores_df, meta_df, by = "SampleID")
            scores_df$LibrarySize <- sample_sums(phyloseq_obj)[scores_df$SampleID]
            
            #Code just to correlate the first axis with the library size
            #rho <- as.numeric(suppressWarnings(cor.test(scores_df$LibrarySize, scores_df$Axis1, method = "spearman")$estimate))
            
            # Determine number of axes to correlate with LibrarySize
            num_axes_available <- sum(grepl("^Axis\\d+$", colnames(scores_df)))
            num_axes_to_correlate <- if (ord_method_upper == "NMDS") {
              min(2, num_axes_available)
            } else {
              min(No_axes_to_correlated_with_LibSize, num_axes_available) #Just make sure you are only uses the avaliable number of axes in case you specified more than are avaliable
            }
            
            # Compute Spearman rho for each axis
            rho_vals <- map_dbl(1:num_axes_to_correlate, function(i) {
              axis_col <- paste0("Axis", i)
              suppressWarnings(cor.test(scores_df$LibrarySize, scores_df[[axis_col]], method = "spearman")$estimate)
            })
            
            # Name the rho values
            names(rho_vals) <- paste0("axis", 1:num_axes_to_correlate)
            
            #Find the variation explained by the first two axes
            var1 <- var2 <- NA
            if (!is.null(var_expl)) {
              var1 <- var_expl[1] * 100
              var2 <- var_expl[2] * 100
            }
            
            #TYPE-based silhouette
            pca_mat <- scores_df[, paste0("Axis", 1:n_components_for_variance)] #Selects the scores up to the number of axes needed to capture the % of variance specified at input
            
            #Calculate silhouette score based on Type only
            scores_df$Type <- factor(scores_df$Type)
            cluster_labels <- as.integer(scores_df$Type)
            sil <- silhouette(cluster_labels, dist(pca_mat))
            sil_df <- as.data.frame(sil[, 1:3])
            sil_df$Type <- scores_df$Type
            colnames(sil_df) <- c("Cluster", "Neighbor", "Silhouette", "Type")
            
            #Find the media sil score per type
            sil_medians <- sil_df %>%
              group_by(Type) %>%
              dplyr::summarise(median_silhouette = median(Silhouette, na.rm = TRUE)) %>%
              pivot_wider(names_from = Type, values_from = median_silhouette, names_prefix = "sil_median_")
            
            #Return results!
            row_result <- c(
              rank_transformation = rank_trans,
              trans_type = trans_type,
              dist_cal_type = dist_type,
              ord_calc_method = ord_method,
              spearman_rho = rho_vals,
              var_expl_axis1 = var1,
              var_expl_axis2 = var2,
              nmds_stress = stress,
              n_components_for_variance = n_components_for_variance  # Add to result
            )
            
            #Save the results in the 
            row_result <- c(row_result, sil_medians)
            results_list[[counter]] <- row_result
            counter <- counter + 1
          })
        }
      }
    }
  }
  
  results_df <- results_list %>%
    map(~{
      .x[setdiff(names(.x), NULL)]  # Drop NULLs
    }) %>%
    bind_rows() %>%
    mutate(across(-c(rank_transformation, trans_type, dist_cal_type, ord_calc_method), as.numeric))
  
  saveRDS(results_df, result_file)
  
  return(list(results_df=results_df))
}

## ---- Run Ordination ----

#Comprehensive ordination plot
run_ordination_with_validation <- function(phyloseq_obj,
                                           rank_transformation = "Genus",
                                           trans_type = "clr",
                                           ord_calc_method = "PCoA",
                                           dist_cal_type = "bray",
                                           component_num = 2,
                                           variable = "Type",
                                           Axis_1 = "Axis1",
                                           Axis_2 = "Axis2",
                                           colors_list = NULL) {

  # Validator
  validate_transformation_distance <- function(trans_type, dist_type, ord_method) {
    ord_method <- toupper(ord_method)
    trans_type <- tolower(trans_type)
    dist_type <- tolower(dist_type)
    
    # Warnings
    
    # CLR requires Euclidean
    if (trans_type == "clr" && dist_type != "euclidean") {
      stop("CLR transformation is only compatible with Euclidean distance.")
    }
    
    # NMDS is discouraged with Euclidean (but allowed)
    if (ord_method == "NMDS" && dist_type == "euclidean") {
      warning("NMDS with Euclidean distance is discouraged in ecological/microbiome data. Prefer Bray or Jaccard.")
    }
    
    # NMDS with CLR is not supported (theoretically questionable)
    if (ord_method == "NMDS" && trans_type == "clr") {
      stop("NMDS is not compatible with CLR-transformed data.")
    }
    
    # Identity (raw counts) with Bray is discouraged
    if (trans_type == "identity" && dist_type == "bray") {
      warning("Using raw counts with Bray-Curtis distance may give misleading results. Consider using relative abundance or log transformation.")
    }
    
    # Identity with Jaccard doesn't make sense
    if (trans_type == "identity" && dist_type == "jaccard") {
      warning("Jaccard distance is intended for presence/absence data. Consider binarizing or using a different transformation.")
    }
    
    # Log transformation with Jaccard is questionable
    if (trans_type == "log" && dist_type == "jaccard") {
      warning("Jaccard distance with log-transformed data is rarely appropriate. Use presence/absence or relative abundance instead.")
    }
    
    # Log transformation with Bray: acceptable, but might want a warning
    if (trans_type == "log" && dist_type == "bray") {
      warning("Log-transformed data with Bray-Curtis distance is valid, but ensure zero handling is appropriate.")
    }
    
    return(dist_type)
  }
  
  # CLR transformation
  clr_transform <- function(mat, pseudocount = 1) {
    mat <- mat + pseudocount
    log_mat <- log(mat)
    gm <- rowMeans(log_mat)
    clr_mat <- sweep(log_mat, 1, gm)
    return(clr_mat)
  }
  
  # Agglomerate
  phylo_trans <- tax_glom(phyloseq_obj, taxrank = rank_transformation)
  otu_mat <- otu_table(phylo_trans)
  if (taxa_are_rows(phylo_trans)) {
    otu_mat <- t(otu_mat)
  }
  
  # Transform
  if (trans_type == "identity") {
    otu_trans <- otu_mat
  } else if (trans_type == "log") {
    otu_trans <- log1p(otu_mat)
  } else if (trans_type == "relab") {
    otu_trans <- sweep(otu_mat, 1, rowSums(otu_mat), FUN = "/")
  } else if (trans_type == "clr") {
    otu_trans <- clr_transform(otu_mat, pseudocount = 1)
  } else {
    stop("Unsupported transformation type.")
  }
  
  # Remove zero-variance columns
  otu_trans <- otu_trans[, apply(otu_trans, 2, var) != 0]
  
  ord_calc_method <- toupper(ord_calc_method)
  dist_type_valid <- validate_transformation_distance(trans_type, dist_cal_type, ord_calc_method)
  
  # Ordination
  if (ord_calc_method == "PCA") {
    dist_type_valid <- "euclidean"
    ord_obj <- prcomp(otu_trans, center = TRUE, scale. = TRUE)
    if (component_num > ncol(ord_obj$x)) {
      stop(paste0("Requested component_num = ", component_num,
                  " but only ", ncol(ord_obj$x), " components are available from PCA."))
    }
    
    scores_df <- ord_obj$x[, 1:component_num] %>% as.data.frame()
    colnames(scores_df) <- paste0("Axis", seq_len(component_num))
    
    var_explained <- ord_obj$sdev^2 / sum(ord_obj$sdev^2)
    
  } else if (ord_calc_method == "PCOA") {
    
    #Binarize for Jaccard
    if (dist_type_valid == "jaccard") {
      dist_mat <- vegdist(otu_trans, method = "jaccard", binary=TRUE)
    } else if (dist_type_valid != "jaccard") {
      dist_mat <- vegdist(otu_trans, method = dist_type_valid)
    }
    
    n_samples <- attr(dist_mat, "Size") 
    ord_obj <- cmdscale(dist_mat, k = n_samples - 1, eig = TRUE)
    scores_df <- as.data.frame(ord_obj$points)
    colnames(scores_df) <- paste0("Axis", 1:ncol(scores_df))
    var_expl <- ord_obj$eig / sum(ord_obj$eig)
    
    
  } else if (ord_calc_method == "NMDS") {
    
    #Binarize for Jaccard
    if (dist_type_valid == "jaccard") {
      dist_mat <- vegdist(otu_trans, method = "jaccard", binary=TRUE)
    } else if (dist_type_valid != "jaccard") {
      dist_mat <- vegdist(otu_trans, method = dist_type_valid)
    }
    
    ord_obj <- metaMDS(dist_mat, trymax = 100)
    scores_df <- as.data.frame(ord_obj$points)
    colnames(scores_df) <- paste0("Axis", 1:ncol(scores_df))
    stress <- ord_obj$stress
  }
  
  # Merge metadata
  scores_df$SampleID <- rownames(scores_df)
  
  metadata_df <- sample_data(phylo_trans) %>% as.matrix() %>% as.data.frame() %>% select(c("SampleID", "Type"))
  
  scores_df <- dplyr::left_join(scores_df, metadata_df, by = "SampleID")


  # Plot
  
  # Label axes with variance explained if available
  x_lab <- Axis_1
  y_lab <- Axis_2
  
  if (exists("var_explained")) {
    axis_nums <- as.numeric(gsub("\\D", "", c(Axis_1, Axis_2)))
    x_lab <- paste0(Axis_1, " (", round(var_explained[axis_nums[1]] * 100, 1), "%)")
    y_lab <- paste0(Axis_2, " (", round(var_explained[axis_nums[2]] * 100, 1), "%)")
  }
  
  p <- ggplot(scores_df, aes(x = !!sym(Axis_1), y = !!sym(Axis_2))) +
    geom_point(aes(fill = !!sym(variable), color = !!sym(variable)), alpha = 1, size = 2, shape = 21) +
    ggtitle(paste0(rank_transformation, " - ", ord_calc_method, " (", trans_type, ", ", dist_cal_type, ")")) +
    stat_ellipse(aes(color = !!sym(variable))) +
    theme_bw() +
    labs(x = x_lab, y = y_lab)+
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = .5),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.text = element_text(size = 12, vjust = 0.5, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(colors_list)) {
    p <- p +
      scale_fill_manual(values = colors_list) +
      scale_color_manual(values = colors_list)
  }
  
  # Axis1 vs library size
  scores_df$LibrarySize <- sample_sums(phyloseq_obj)[scores_df$SampleID]
  lib_size_with_Axis1 <- ggplot(scores_df, aes(x = LibrarySize, y = Axis1)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(x = "Library Size", y = paste0(Axis_1), title = paste0(Axis_1, " vs Library Size")) +
    theme_minimal()
  lib_size_with_Axis1_corr <- cor.test(scores_df$LibrarySize, scores_df[[Axis_1]], method = "spearman")
  
  
  # Scree/L-plot for PCA or PCoA
  if (ord_calc_method %in% c("PCA", "PCOA")) {
    if (ord_calc_method == "PCA") {
      eig_vals <- ord_obj$sdev^2
    } else if (ord_calc_method == "PCOA") {
      eig_vals <- ord_obj$eig
    }
    
    var_explained <- eig_vals / sum(eig_vals)
    axis_labels <- paste0("Axis", seq_along(var_explained))
    scree_df <- data.frame(Axis = factor(axis_labels, levels = axis_labels),
                           VarianceExplained = var_explained)
    
    l_plot <- ggplot(scree_df[1:component_num, ], aes(x = Axis, y = VarianceExplained)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = round(VarianceExplained, 3)), vjust = -0.5, size = 3) +
      labs(title = paste("Variance Explained per Axis -", ord_calc_method),
           y = "Proportion of Variance Explained", x = NULL) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    l_plot <- NULL
  }
  
  return(list(
    plot = p,
    scores = scores_df,
    ordination = ord_obj,
    lib_size_with_Axis1_plot = lib_size_with_Axis1, 
    lib_size_with_Axis1_corr = lib_size_with_Axis1_corr,
    l_plot = l_plot
  ))
}

## ---- Analyze Silhouette Score ----

# Analyze Silhouette score using "Type" as cluster label and return plot
analyze_type_clustering_on_pca <- function(scores_df,
                                           component_num = 2,
                                           colors_list = NULL) {

  # Extract PCA axes
  pca_matrix <- scores_df[, paste0("Axis", 1:component_num)]
  
  # Ensure 'Type' is a factor
  scores_df$Type <- factor(scores_df$Type)
  
  # Calculate distance matrix and silhouette scores using "Type" as cluster labels
  dist_matrix <- dist(pca_matrix)
  cluster_labels <- as.integer(scores_df$Type)  # silhouette() needs numeric cluster labels
  
  sil <- silhouette(cluster_labels, dist_matrix)
  
  sil_df <- as.data.frame(sil[, 1:3])
  sil_df$Type <- scores_df$Type
  colnames(sil_df) <- c("Cluster", "Neighbor", "Silhouette", "Type")
  
  # Boxplot of silhouette scores by Type
  sil_plot <- ggplot(sil_df, aes(x = Type, y = Silhouette, fill = Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 1) +
    geom_jitter(width = 0.2, size = 1, alpha = 1) +
    stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                 vjust = .5, color = "white", fontface = "bold", size = 4.5) +
    labs(title = "Silhouette Scores by Type",
         x = "Type",
         y = "Silhouette Score") +
    theme_bw() +
    theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if (!is.null(colors_list)) {
    sil_plot <- sil_plot + scale_fill_manual(values = colors_list)
  }
  
  return(list(
    scores = scores_df,
    silhouette_df = sil_df,
    sil_plot = sil_plot
  ))
}



# ---- DA Analysis ----

## ---- Maaslin2 ----

### Iterate Maaslin2 Params
#Example:
#resolutions <- c("Species", "Genus", "Phylum")
#iterate_maaslin2(ps_obj = ps_obj_ped_fs, 
#                 iterative_methods = iterative_methods,
#                 resolutions = resolutions,
#                 group = "condition",
#                 qval_threshold = .25, 
#                 plot_colors = maaslin2_colors,
#                 percentage=TRUE)

iterate_maaslin2 <- function(ps_obj, iterative_methods, resolutions, group, qval_threshold, percentage, plot_colors) {
  
  # Define the options for each parameter
  analysis_methods <- c("LM", "CPLM", "ZINB", "NEGBIN")
  transforms <- c("NONE", "LOG", "LOGIT", "AST")
  normalizations <- c("TSS", "TMM", "CSS", "CLR")
  enriched_taxa <- list()
  
  # Generate all combinations using expand.grid
  combinations <- expand.grid(analysis_method = analysis_methods,
                              transform = transforms,
                              normalization = normalizations)
  
  # Convert each row into a list and store in a list
  iterative_methods <- apply(combinations, 1, function(x) {
    list(analysis_method = x["analysis_method"],
         transform = x["transform"],
         normalization = x["normalization"])
  })
  
  # Create a unique hash for the parameters
  param_hash <- digest(list(ps_obj, resolutions, group, qval_threshold))
  # File path for the analysis results
  result_file <- file.path("..saved_analysis_files/", paste0("iterative_maaslin2_result_", param_hash, ".rds"))
  
  if (file.exists(result_file)) {
    message("Analysis already run. Loading results...")
    summary_table <- readRDS(result_file)
  } else {
    message("Running analysis...")
    
    # Initialize an empty data frame to store results
    summary_table <- data.frame(
      analysis_method = character(),
      transform = character(),
      normalization = character(),
      resolution = character(),
      significant_features = integer(),
      stringsAsFactors = FALSE
    )
    
    # Loop through each resolution (taxa_level)
    for (taxa_level in resolutions) {
      print(taxa_level)
      
      # Prepare counts and metadata input
      counts_input <- ps_obj %>% 
        tax_glom(taxa_level) %>% 
        psmelt() %>%
        dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
        pivot_wider(names_from = taxa_level, values_from = "Abundance") %>%
        dplyr::select(-group) %>% 
        column_to_rownames(var = "Sample")
      
      metadata_input <- sample_data(ps_obj) %>% as.matrix() %>% as.data.frame()
      
      # Loop through each combination of parameters
      for (params in iterative_methods) {
        
        # Extract parameters for the current iteration
        analysis_method <- params$analysis_method
        transform <- params$transform
        normalization <- params$normalization
        
        # Run Maaslin2 with current parameters inside tryCatch
        tryCatch({
          maaslin2_results <- Maaslin2(
            counts_input,
            metadata_input,
            fixed_effects = group,
            plot_heatmap = FALSE,
            plot_scatter = FALSE,
            output = "../saved_analysis_files/maaslin_iterations",
            analysis_method = analysis_method,
            normalization = normalization,
            transform = transform
          )
          
          # Filter results based on q-value threshold
          maaslin_res_filt <- maaslin2_results$results %>% 
            as.data.frame() %>%
            filter(qval <= qval_threshold)
          
          # Count unique significant features
          significant_count <- length(unique(maaslin_res_filt$feature))
          features_list <- unique(maaslin_res_filt$feature)
          features_list_str <- paste(features_list, collapse = ", ")
          
          
          # Append results to the summary table
          summary_table <- rbind(summary_table, data.frame(
            analysis_method = analysis_method,
            transform = transform,
            normalization = normalization,
            resolution = taxa_level,
            significant_features = significant_count,
            features_list = features_list_str,
            stringsAsFactors = FALSE
          ))
          
        }, error = function(e) {
          # Print an error message and skip to the next iteration
          message(paste("Error in Maaslin2 for parameters:",
                        "analysis_method =", analysis_method,
                        "transform =", transform,
                        "normalization =", normalization))
          message("Error details:", e$message)
        })
        
        
      }
    }
    # Save the results after each iteration
    saveRDS(summary_table, result_file)
  }
  
  # Return the summary table
  summary_table_features_list <- summary_table %>%
    mutate(method = paste(resolution, analysis_method, normalization, transform, sep = "_")) %>%
    select(c("resolution", "method", "features_list"))
  
  summary_table <- summary_table %>% select(-c("features_list"))
  
  summary_table$resolution = factor(summary_table$resolution, 
                                    levels = resolutions)
  
  
  if (percentage == TRUE) {
    
    #Get the number of unique taxa for the proportion calc
    total_sp <- length(unique(tax_table(ps_obj)[, "Species"]))
    print(paste0("Total species: ", total_sp))
    total_gn <- length(unique(tax_table(ps_obj)[, "Genus"]))
    print(paste0("Total genera: ", total_gn))
    total_phy <- length(unique(tax_table(ps_obj)[, "Phylum"]))
    print(paste0("Total phyla: ", total_phy))
    
    summary_table  <- summary_table %>%
      mutate(normalization_transform = paste(normalization, transform, sep = "_")) %>%
      mutate(percent_significant_features = case_when(
        resolution == "Phylum" ~ round(((significant_features/total_phy)*100), 0),
        resolution == "Genus" ~ round(((significant_features/total_gn)*100), 0),
        resolution == "Species" ~ round(((significant_features/total_sp)*100), 0))) 
    
    p <- summary_table %>%
      ggplot(., aes(x = analysis_method, y = percent_significant_features, fill = normalization_transform)) +
      geom_col(position = position_dodge(width = .8), width = 0.5, color="black") +
      geom_text(aes(label = percent_significant_features), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, size = 2) +  # Ensure labels match bar positions      theme_bw(base_size = 14) +
      facet_wrap(~resolution, scales = "free", ncol=1) +
      theme_bw(base_size = 14) +
      theme(strip.text.x = element_text(size = 16))+
      labs(
        title = "Significant Features by Analysis Method, Resolution, and Normalization",
        x = "Analysis Method",
        y = "% of Total Features Found Signficiant"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_fill_manual(values = plot_colors) + 
      ylim(0, 100) 
    
  }else if (percentage == FALSE) {
    
    p <- summary_table %>%
      mutate(normalization_transform = paste(normalization, transform, sep = "_")) %>%
      ggplot(., aes(x = analysis_method, y = significant_features, fill = normalization_transform)) +
      geom_col(position = position_dodge(width = .8), width = 0.5, color="black") +
      geom_text(aes(label = significant_features), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, size = 2) +  # Ensure labels match bar positions
      facet_wrap(~resolution, scales = "free", ncol=1) +
      theme(strip.text.x = element_text(size = 16))+
      theme_bw(base_size = 14) +
      labs(
        title = "Significant Features by Analysis Method, Resolution, and Normalization",
        x = "Analysis Method",
        y = "Significant Features"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_fill_manual(values = plot_colors)
  }
  return(list(plot = p, results = summary_table))  
}


# Function to do the same iterrations, but with the picrust2 output
iterate_maaslin2_picrust2 <- function(counts_input, metadata, iterative_methods, group, qval_threshold, percentage, plot_colors) {
  
  # Define the options for each parameter
  analysis_methods <- c("LM", "CPLM", "ZINB", "NEGBIN")
  transforms <- c("NONE", "LOG", "LOGIT", "AST")
  normalizations <- c("TSS", "TMM", "CSS", "CLR")
  enriched_taxa <- list()
  
  # Generate all combinations using expand.grid
  combinations <- expand.grid(analysis_method = analysis_methods,
                              transform = transforms,
                              normalization = normalizations)
  
  # Convert each row into a list and store in a list
  iterative_methods <- apply(combinations, 1, function(x) {
    list(analysis_method = x["analysis_method"],
         transform = x["transform"],
         normalization = x["normalization"])
  })
  
  # Create a unique hash for the parameters
  param_hash <- digest(list(counts_input, group, qval_threshold))
  # File path for the analysis results
  result_file <- file.path("..saved_analysis_files/", paste0("iterative_maaslin2_picrust2_result_", param_hash, ".rds"))
  
  
  # Initialize an empty data frame to store results
  summary_table <- data.frame(
    analysis_method = character(),
    transform = character(),
    normalization = character(),
    significant_features = integer(),
    stringsAsFactors = FALSE
  )
  
  
  # Loop through each combination of parameters
  for (params in iterative_methods) {
    
    # Extract parameters for the current iteration
    analysis_method <- params$analysis_method
    transform <- params$transform
    normalization <- params$normalization
    
    # Run Maaslin2 with current parameters inside tryCatch
    tryCatch({
      maaslin2_results <- Maaslin2(
        counts_input,
        metadata,
        fixed_effects = group,
        plot_heatmap = FALSE,
        plot_scatter = FALSE,
        output = "../saved_analysis_files/maaslin_iterations",
        analysis_method = analysis_method,
        normalization = normalization,
        transform = transform
      )
      
      # Filter results based on q-value threshold
      maaslin_res_filt <- maaslin2_results$results %>% 
        as.data.frame() %>%
        filter(qval <= qval_threshold)
      
      # Count unique significant features
      significant_count <- length(unique(maaslin_res_filt$feature))
      features_list <- unique(maaslin_res_filt$feature)
      features_list_str <- paste(features_list, collapse = ", ")
      
      
      # Append results to the summary table
      summary_table <- rbind(summary_table, data.frame(
        analysis_method = analysis_method,
        transform = transform,
        normalization = normalization,
        significant_features = significant_count,
        features_list = features_list_str,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      # Print an error message and skip to the next iteration
      message(paste("Error in Maaslin2 for parameters:",
                    "analysis_method =", analysis_method,
                    "transform =", transform,
                    "normalization =", normalization))
      message("Error details:", e$message)
    })
    
    
  }
  
  
  if (percentage == TRUE) {
    
    #Get the number of unique taxa for the proportion calc
    total_path <- nrow(counts_input)
    print(paste0("Total pathways: ", total_path))
    
    summary_table  <- summary_table %>%
      mutate(normalization_transform = paste(normalization, transform, sep = "_")) %>%
      mutate(percent_significant_features = round(((significant_features/total_path)*100), 0)) 
    
    p <- summary_table %>%
      ggplot(., aes(x = analysis_method, y = percent_significant_features, fill = normalization_transform)) +
      geom_col(position = position_dodge(width = .8), width = 0.5, color="black") +
      geom_text(aes(label = percent_significant_features), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, size = 2) +  # Ensure labels match bar positions      theme_bw(base_size = 14) +
      theme_bw(base_size = 14) +
      theme(strip.text.x = element_text(size = 16))+
      labs(
        title = "Significant Features by Analysis Method, Resolution, and Normalization",
        x = "Analysis Method",
        y = "% of Total Features Found Signficiant"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_fill_manual(values = plot_colors) + 
      ylim(0, 100) 
    
  }else if (percentage == FALSE) {
    
    p <- summary_table %>%
      mutate(normalization_transform = paste(normalization, transform, sep = "_")) %>%
      ggplot(., aes(x = analysis_method, y = significant_features, fill = normalization_transform)) +
      geom_col(position = position_dodge(width = .8), width = 0.5, color="black") +
      geom_text(aes(label = significant_features), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, size = 2) +  # Ensure labels match bar positions
      theme(strip.text.x = element_text(size = 16))+
      theme_bw(base_size = 14) +
      labs(
        title = "Significant Features by Analysis Method, Resolution, and Normalization",
        x = "Analysis Method",
        y = "Significant Features"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_fill_manual(values = plot_colors)
  }
  return(list(plot = p, results = summary_table))  
}


### Run Maaslin2
#Example:
#run_Maaslin2(ps_obj = ps_obj_ped_fs, 
#                                 taxa_level = "Genus", 
#                                 group = "condition",
#                                 analysis_method = "NEGBIN", 
#                                 normalization = "TMM",
#                                 transform = "NONE", 
#                                 plot_colors = plot_colors,
#                                 plot_type = "box", 
#                                 qval_threshold= 0.1, 
#                                 plot_heights = c(2, 9))

run_Maaslin2 <- function(ps_obj, taxa_level, group, analysis_method, normalization, transform, plot_colors, plot_type, qval_threshold){
  
  # Create a unique hash for the parameters
  param_hash <- digest(list(ps_obj, taxa_level, group, analysis_method, normalization, transform, plot_colors, plot_type, qval_threshold))
  # File path for the analysis results
  result_file <- file.path("../saved_analysis_files/", paste0("maaslin2_result_", param_hash, ".rds"))
  
  if (file.exists(result_file)) {
    message("Analysis already run. Loading results...")
    maaslin2_results <- readRDS(result_file)
  } else {
    message("Running analysis...")
    
    counts_input <- ps_obj %>% tax_glom(., taxa_level) %>% psmelt() %>%
      dplyr::select(c("Sample", "Abundance",taxa_level, group)) %>%
      pivot_wider(names_from=taxa_level, values_from="Abundance") %>%
      dplyr::select(-c(group)) %>% 
      column_to_rownames(var="Sample") 
    
    metadata_input <- sample_data(ps_obj) %>% as.matrix() %>% as.data.frame()
    
    maaslin2_results <- tryCatch({
      Maaslin2(
        input_data = counts_input,
        input_metadata = metadata_input,
        fixed_effects = group,
        plot_heatmap = FALSE,
        plot_scatter = FALSE,
        output = paste0("../saved_analysis_files/", Sys.Date(), "_maaslin2_output"),
        analysis_method = analysis_method,
        normalization = normalization,
        transform = transform
      )
    }, error = function(e) {
      message("Maaslin2 failed: ", e$message)
      return(NULL)
    })
    
    # If run failed
    if (is.null(maaslin2_results)) {
      return(invisible(NULL))
    }
    
    # Save the results
    saveRDS(maaslin2_results, result_file)
  }
  
 
  #PLOT EVERYTHING
  
  maaslin_res_filt <- maaslin2_results$results %>% as.data.frame() %>% filter(qval <= qval_threshold) %>%
    dplyr::mutate(Enrichment = ifelse(coef > 0, "Plaque", "Abscess"))
  
  # Assign the results df to the global environment to output
  assign(paste0("maaslin2_results_", taxa_level, "_", qval_threshold), maaslin_res_filt, envir = .GlobalEnv)
  
  if (nrow(maaslin_res_filt) == 0) {
    stop("There are no results to display.")
  }
  
  df_fig <- maaslin_res_filt %>%
    mutate(feature = gsub("g__", "", feature)) %>%
    mutate(cols=ifelse(coef > 0, "CNTRL_enriched", "ABNORM_enriched")) %>%
    arrange(Enrichment, abs(coef)) %>% 
    mutate(feature_order = factor(feature, levels = unique(feature)))
  
  #Barplot
  bar_plot <- df_fig %>%
    ggplot(aes(x = feature_order, y = abs(coef), fill = Enrichment)) +
    geom_bar(stat = "identity", color="black") +
    scale_fill_manual(values = plot_colors) +
    coord_flip() +
    labs(title = "Maaslin2 Results", x=NA, y = "coef") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          axis.text.y = element_text(face="italic"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  dot_plot <- df_fig %>%
    ggplot(aes(x = abs(coef), y = feature_order)) +
    geom_point(aes(size = -log10(pval), color = Enrichment, fill = Enrichment)) + # Map size and color
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_size(range = c(1, 5), name = "-log10(P-value)") + # Adjust dot size and legend label
    labs(title = NULL, x="Log fold change", y = NULL, color = "Enrichment") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
      plot.title = element_text(hjust = 0.5),  # Center the title
      axis.text.y = element_text(hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  
  #CREATE PLOTS OF COUNTS, TSS, and CLR
  physeq_glom_df_counts <- tax_glom(ps_obj, taxa_level) %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
    mutate(Count_Type = "Count")
  
  #Create Relative Abundance Table
  physeq_glom_df_relab <- tax_glom(ps_obj, taxa_level) %>%
    transform_sample_counts(., function(x) (x + 0.000001 - min(x))) %>%
    transform_sample_counts(., function(x) (x) / sum(x)) %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
    mutate(Count_Type = "RelAb")
  
  #Create CLR table
  physeq_glom_df_clr<- tax_glom(ps_obj, taxa_level) %>%
    microbiome::transform(., transform = "clr", target="sample") %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level, group)) %>%
    mutate(Count_Type = "CLR")    
  
  merged <- rbind(physeq_glom_df_counts, physeq_glom_df_relab, physeq_glom_df_clr)
  
  # Filter for taxa in Maaslin2 results
  filtered_data <- merged %>%
    dplyr::filter(!!sym(taxa_level) %in% maaslin_res_filt$feature) %>%
    mutate(Count_Type = factor(Count_Type, levels = c("Count", "RelAb", "CLR"))) # Reorder levels
  
  # Create dot plot
  if (plot_type == "dot") {
    plot <- filtered_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_point(position = position_jitter(seed = 1, width = 0.2)) +
      facet_wrap(vars(!!sym(taxa_level), Count_Type), ncol=3, scales = "free") + # Use facet_wrap with free scales
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    #Box plot
  }else if (plot_type == "box") {
    plot <- filtered_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_boxplot(alpha = 0.5, outlier.shape = 8) + # Add boxplot
      #geom_jitter(width = 0.2, alpha = 0.7, size = 1.5) + # Add jittered points
      facet_wrap(vars(!!sym(taxa_level), Count_Type), ncol = 3, scales = "free") + # Facet
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    #Violin
  }else if (plot_type == "violin") {
    plot <- filtered_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_violin(alpha = 0.5, scale = "width") + # Add violin plot with partial transparency
      geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 1.5, alpha = 0.7) + # Add jittered points
      facet_wrap(vars(!!sym(taxa_level), Count_Type), ncol = 3, scales = "free") + # Facet
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    stop("Invalid plot type. Please specify 'box', 'violin', 'density' or 'dot'.")
  }
  
  
  return(list(results = df_fig, bar_plot = bar_plot, dot_plot = dot_plot))
  
}

## ---- ANCOM ----

### Run ANCOM
#Example:

#runANCOM(merged_phylo_obj_fs, 
#         tax_level = "Genus", 
#         group = "condition", 
#         Log2FC_cutoff = 1, 
#         name_of_saved_results = "ancombc_results_genus", 
#         plot_type = "dot", 
#         plot_heights = c(1, 8))

runANCOM <- function(ps_obj, tax_level, group, Log2FC_cutoff=NULL, name_of_saved_results=NULL, plot_type, plot_colors){
  
  # Create a unique hash for the parameters
  param_hash <- digest(list(ps_obj, tax_level, group, Log2FC_cutoff=NULL))
  # File path for the analysis results
  result_file <- file.path("../saved_analysis_files/", paste0("ancombc2_result_", param_hash, ".rds"))
  
  if (file.exists(result_file)) {
    message("Analysis already run. Loading results...")
    ancombc2_results <- readRDS(result_file)
  } else {
    message("Running analysis...")
    
    # Run ANCOMBC2
    ancombc2_results <- tryCatch({
      res <- ancombc2(
        data = ps_obj,
        assay_name = "counts",
        tax_level = tax_level,
        fix_formula = group,
        p_adj_method = "holm",
        group = group,
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        n_cl = 1
      )
      saveRDS(res, result_file)
      res
    }, error = function(e) {
      message("ANCOMBC2 failed: ", e$message)
      return(NULL)
    })
  
  }
  df_fig <- ancombc2_results$res
  
  
  #Save the column names to reference below
  #This seems really complicated but it's not, it just makes the function flexible to whatever group you are comparing.
  
  #DFF column
  dff_col_name <- grep(paste0("^diff_", group), colnames(df_fig), value = TRUE)
  # If you want to check and store it
  if (length(dff_col_name) > 0) {
    # Save the first matching column name (if multiple matches)
    dff_col_name <- dff_col_name[1]
  } else {
    print("No matching column found.")
  }
  
  #LFC column
  lfc_col_name <- grep(paste0("^lfc_", group), colnames(df_fig), value = TRUE)
  # If you want to check and store it
  if (length(lfc_col_name) > 0) {
    # Save the first matching column name (if multiple matches)
    lfc_col_name <- lfc_col_name[1]
  } else {
    print("No matching column found.")
  }
  
  #SS condition column
  ss_col_name <- grep(paste0("^passed_ss_", group), colnames(df_fig), value = TRUE)
  # If you want to check and store it
  if (length(ss_col_name) > 0) {
    # Save the first matching column name (if multiple matches)
    ss_col_name <- ss_col_name[1]
  } else {
    print("No matching column found.")
  }
  
  #SE condition column
  se_col_name <- grep(paste0("^se_", group), colnames(df_fig), value = TRUE)
  # If you want to check and store it
  if (length(se_col_name) > 0) {
    # Save the first matching column name (if multiple matches)
    se_col_name <- se_col_name[1]
  } else {
    print("No matching column found.")
  }
  
  #P-val condition column
  p_col_name <- grep(paste0("^p_", group), colnames(df_fig), value = TRUE)
  # If you want to check and store it
  if (length(p_col_name) > 0) {
    # Save the first matching column name (if multiple matches)
    p_col_name <- p_col_name[1]
  } else {
    print("No matching column found.")
  }
  
  #Create a dataframe for results
  df_fig <- ancombc2_results$res %>%
    dplyr::filter(!!sym(dff_col_name) == TRUE) %>% 
    dplyr::arrange(desc(!!sym(lfc_col_name))) %>%
    dplyr::mutate(Enrichment = ifelse(!!sym(lfc_col_name) > 0, "Plaque", "Abscess"),
                  color = ifelse(!!sym(ss_col_name) == 1, "aquamarine3", "black"))  # Color text based on ss condition
  
  df_fig$taxon = factor(df_fig$taxon, levels = df_fig$taxon)
  df_fig$Enrichment = factor(df_fig$Enrichment, 
                             levels = c("Plaque", "Abscess"))
  
  # Assign the results df to the global environment to output
  assign(paste0("ancombc2_results_df_", tax_level),df_fig, envir = .GlobalEnv)
  
  # Check if the dataframe is empty
  if (nrow(df_fig) == 0) {
    message("ANCOMBC2 ran successfully but found no significant taxa.")
    return(invisible(NULL))  # Don't error, just end the function quietly
  }
  

  #DOT PLOT
  # Conditionally filter the data based on Log2FC_cutoff
  if (!is.null(Log2FC_cutoff)) {
    df_fig <- df_fig %>% filter(abs(!!sym(lfc_col_name)) > Log2FC_cutoff)
  }
  
  dot_plot <- df_fig %>%
    ggplot(aes(x = abs(!!sym(lfc_col_name)), y = taxon)) +
    geom_point(aes(size = -log10(!!sym(p_col_name)), color = Enrichment)) + # Map size and color
    scale_color_manual(values = plot_colors) +
    scale_size(range = c(1, 5), name = "-log10(P-value)") + # Adjust dot size and legend label
    labs(title = NULL, x="Log fold change", y = NULL, color = "Enrichment") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
      plot.title = element_text(hjust = 0.5),  # Center the title
      axis.text.y = element_text(hjust = 1, color = df_fig$color),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  
  #Bar plot
  df_fig$taxon = factor(df_fig$taxon, levels = df_fig$taxon)
  df_fig$Enrichment = factor(df_fig$Enrichment, 
                             levels = c("Plaque", "Abscess"))
  
  bar_plot <- df_fig %>%
    ggplot(aes(x = taxon, y = abs(!!sym(lfc_col_name)), fill = Enrichment)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = abs(!!sym(lfc_col_name)) - abs(!!sym(se_col_name)), ymax = abs(!!sym(lfc_col_name)) + abs(!!sym(se_col_name))), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = NULL) + 
    scale_fill_manual(values = plot_colors) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_text(hjust = 1, color = df_fig$color),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_flip()
  
  
  if (!is.null(name_of_saved_results)) {
    assign(name_of_saved_results, df_fig, envir = .GlobalEnv)
  }
  
  
  #CREATE PLOTS OF COUNTS, TSS, and CLR
  physeq_glom_df_counts <- tax_glom(ps_obj, tax_level) %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", tax_level, group)) %>%
    mutate(Count_Type = "Count")
  
  #Create Relative Abundance Table
  physeq_glom_df_relab <- tax_glom(ps_obj, tax_level) %>%
    transform_sample_counts(., function(x) (x + 0.000001 - min(x))) %>%
    transform_sample_counts(., function(x) (x) / sum(x)) %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", tax_level, group)) %>%
    mutate(Count_Type = "RelAb")
  
  
  #Create CLR table
  physeq_glom_df_clr<- tax_glom(ps_obj, tax_level) %>%
    microbiome::transform(., transform = "clr", target="sample") %>% psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", tax_level, group)) %>%
    mutate(Count_Type = "CLR")    
  
  merged <- rbind(physeq_glom_df_counts, physeq_glom_df_relab, physeq_glom_df_clr)
  
  # Filter for taxa in ANCOM results
  filtered_data <- merged %>%
    filter(!!sym(tax_level) %in% df_fig$taxon) %>%
    mutate(Count_Type = factor(Count_Type, levels = c("Count", "RelAb", "CLR"))) # Reorder levels
  
  # Create dot plot
  if (plot_type == "dot") {
    plot <- filtered_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_point(position = position_jitter(seed = 1, width = 0.2)) +
      facet_wrap(vars(!!sym(tax_level), Count_Type), ncol=3, scales = "free") + # Use facet_wrap with free scales
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    #Box plot
  }else if (plot_type == "box") {
    plot <- filtered_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_boxplot(alpha = 0.5, outlier.shape = 8) + # Add boxplot
      #geom_jitter(width = 0.2, alpha = 0.7, size = 1.5) + # Add jittered points
      facet_wrap(vars(!!sym(tax_level), Count_Type), ncol = 3, scales = "free") + # Facet
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    #Violin
  }else if (plot_type == "violin") {
    plot <- filtered_data %>%
      ggplot(aes(x = !!sym(group), y = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_violin(alpha = 0.5, scale = "width") + # Add violin plot with partial transparency
      geom_jitter(position = position_jitter(seed = 1, width = 0.2), size = 1.5, alpha = 0.7) + # Add jittered points
      facet_wrap(vars(!!sym(tax_level), Count_Type), ncol = 3, scales = "free") + # Facet
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    #Density
  }else if (plot_type == "density") {
    plot <- filtered_data %>%
      ggplot(aes(x = Abundance, color = !!sym(group), fill = !!sym(group))) + 
      geom_density(alpha = 0.5) + # Add density plot with partial transparency
      facet_wrap(vars(!!sym(tax_level), Count_Type), ncol = 3, scales = "free") + # Facet
      theme_bw(base_size = 10) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    
  } else {
    stop("Invalid plot type. Please specify 'box', 'violin', 'density' or 'dot'.")
  }
  
  return(list(bar_plot = bar_plot, dot_plot = dot_plot, full_results = ancombc2_results$res, filtered_res = df_fig))


}


## ---- Aldex2 ----


### Iterate ALDEX
iterate_aldex2 <- function(ps_obj, iterative_methods, resolutions, group, plot_colors, percentage) {
  
  #p_adjusts = c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY")
  # Define the options for each parameter
  transforms <- c("identity", "log10", "log10p")
  analysis_methods = c("t.test", "wilcox.test", "kruskal", "glm_anova")
  normalization = c("none", "rarefy", "TSS", "TMM", "RLE", "CSS", "CLR", "CPM")
  denoms = c("all", "iqlr", "zero", "lvha")
  paired = c("TRUE", "FALSE")
  
  enriched_taxa <- list()
  
  # Generate all combinations using expand.grid
  combinations <- expand.grid(analysis_method = analysis_methods,
                              transform = transforms,
                              normalization = normalization,
                              denoms = denoms, 
                              paired=paired)

    # Convert each row into a list and store in a list
  iterative_methods <- apply(combinations, 1, function(x) {
    list(analysis_method = x["analysis_method"],
         transform = x["transform"],
         normalization = x["normalization"],
         denoms = x["denoms"],
         paired = x["paired"])})


  # Create a unique hash for the parameters
  param_hash <- digest(list( ps_obj, resolutions, group))
  # File path for the analysis results
  result_file <- file.path("../saved_analysis_files/", paste0("iterative_aldex2_result_", param_hash, ".rds"))

  if (file.exists(result_file)) {
    message("Analysis already run. Loading results...")
    summary_table <- readRDS(result_file)
    
  } else {
    message("Running analysis...")
  
    # Initialize an empty data frame to store results
    summary_table <- data.frame(
      analysis_method = character(),
      transform = character(),
      normalization = character(),
      resolution = character(),
      denom = character(),
      paired = character(),
      significant_features = integer(),
      features_list = character(),
      stringsAsFactors = FALSE
    )
    
    
    # Loop through each resolution (taxa_level)
    for (taxa_level in resolutions) {
      print(taxa_level)
      
    
     # Loop through each combination of parameters
          for (params in iterative_methods) {
            
            # Extract parameters for the current iteration
            analysis_method <- as.character(params$analysis_method)
            transform <- as.character(params$transform)
            normalization <- as.character(params$normalization)
            denom <- as.character(params$denom)
            paired <- as.character(params$paired)
            taxa <- as.character(taxa_level)
            
            print(paste0(analysis_method, transform, normalization, denom, paired, taxa))
            
        tryCatch({
        
            aldex2_results <- microbiomeMarker::run_aldex(
                ps_obj,
                group,
                taxa_rank = taxa_level,
                transform = transform,
                norm = normalization,
                method = analysis_method,
                p_adjust = "BH",
                pvalue_cutoff = 0.05,
                mc_samples = 128,
                denom = denom,
                paired = paired)
            
            aldex2_res <- aldex2_results@marker_table %>% as.data.frame()
            significant_count <- length(unique(aldex2_res$feature))
            features_list <- unique(aldex2_res$feature)
            features_list_str <- paste(features_list, collapse = ", ")
        
            summary_table <- rbind(summary_table, data.frame(
                analysis_method = analysis_method,
                transform = transform,
                normalization = normalization,
                resolution = taxa_level,
                denom = denom, 
                paired = paired,
                significant_features = significant_count,
                features_list = features_list_str,
                stringsAsFactors = FALSE
            ))
        }, error = function(e) {
                message(
                    "An error occurred during the analysis with the following parameters:",
                    paste("  Analysis Method:", analysis_method),
                    paste("  Transform:", transform),
                    paste("  Normalization:", normalization),
                    paste("  Denominator:", denom),
                    paste("  Paired:", paired),
                    paste("  Resolution:", taxa_level),
                    paste("Error Details:", conditionMessage(e)),
                    sep = "\n"
            )
        })

      }
    }
      # Save the results after each iteration
      saveRDS(summary_table, result_file)
  }
  
    #Output
  assign("aldex2_iterative_summary_table", summary_table, envir = .GlobalEnv)
  
   # Return the summary table
  summary_table_features_list <- summary_table %>%
    mutate(method = paste(normalization, transform, denoms, paired,  sep = "_")) %>%
    select(c("resolution", "method", "features_list"))
  
  summary_table <- summary_table %>% select(-c("features_list"))

  summary_table$resolution = factor(summary_table$resolution, 
                               levels = resolutions)
  
  
  if (percentage == TRUE) {
    
  #Get the number of unique taxa for the proportion calc
  total_sp <- length(unique(tax_table(ps_obj)[, "Species"]))
  print(paste0("Total species: ", total_sp))
  total_gn <- length(unique(tax_table(ps_obj)[, "Genus"]))
  print(paste0("Total genera: ", total_gn))
  total_phy <- length(unique(tax_table(ps_obj)[, "Phylum"]))
  print(paste0("Total phyla: ", total_phy))

  summary_table  <- summary_table %>%
    mutate(normalization_transform = paste(normalization, transform, denom, paired, sep = "_")) %>%
    mutate(percent_significant_features = case_when(
                resolution == "Phylum" ~ round(((significant_features/total_phy)*100), 0),
                resolution == "Genus" ~ round(((significant_features/total_gn)*100), 0),
                resolution == "Species" ~ round(((significant_features/total_sp)*100), 0))) %>%
    filter(significant_features > 0)
  
    p <- summary_table %>%
    ggplot(., aes(x = analysis_method, y = percent_significant_features, fill = normalization_transform)) +
      geom_col(position = position_dodge(width = .8), width = 0.5, color="black") +
    geom_text(aes(label = percent_significant_features), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 2) +  # Ensure labels match bar positions      theme_bw(base_size = 14) +
    facet_wrap(~resolution, scales = "free", ncol=1) +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 16))+
    labs(
        title = "Significant Features by Analysis Method, Resolution, and Normalization",
        x = "Analysis Method",
        y = "% of Total Features Found Signficiant"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #scale_fill_manual(values = plot_colors) + 
    #ylim(0, 100) 
  
  }else if (percentage == FALSE) {
  
  summary_table <- summary_table %>%
    mutate(normalization_transform = paste(normalization, transform, denom, paired,  sep = "_")) %>%
    filter(significant_features > 0)
  
  assign("TEST", summary_table, envir = .GlobalEnv)

  p <- summary_table %>% ggplot(., aes(x = analysis_method, y = significant_features, fill = normalization_transform)) +
      geom_col(position = position_dodge(width = .8), width = 0.5, color="black") +
    geom_text(aes(label = significant_features), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 2) +  # Ensure labels match bar positions
      facet_wrap(~resolution, scales = "free", ncol=1) +
     theme(strip.text.x = element_text(size = 16))+
     theme_bw(base_size = 14) +
      labs(
        title = "Significant Features by Analysis Method, Resolution, and Normalization",
        x = "Analysis Method",
        y = "Significant Features"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
   # scale_fill_manual(values = plot_colors)
  }

  return(list(results = summary_table, plot = p))
}



### Run Aldex

run_aldex2 <- function(ps_obj, group, tax_rank, method, transform, normalization, plot_colors) {
  library(dplyr)
  library(ggplot2)
  library(microbiomeMarker)
  library(patchwork)
  
  message("Running ALDEx2...")
  
  # Try running ALDEx2
  aldex2_results <- tryCatch({
    microbiomeMarker::run_aldex(
      ps_obj,
      group,
      method = method,
      norm = normalization,
      transform = transform,
      taxa_rank = tax_rank,
      pvalue_cutoff = 0.05,
      denom = "all",
      paired = FALSE
    )
  }, error = function(e) {
    message("ALDEx2 failed: ", e$message)
    return(NULL)
  })
  
  # Exit if it failed
  if (is.null(aldex2_results)) {
    return(invisible(NULL))
  }
  
  # Extract results
  if (!"marker_table" %in% slotNames(aldex2_results) || nrow(aldex2_results@marker_table) == 0) {
    message("ALDEx2 returned no marker results.")
    return(invisible(NULL))
  }
  
  # Convert and assign result table
  aldex2_res <- aldex2_results@marker_table %>% as.matrix() %>% as.data.frame()
  assign(paste0("aldex2_res_", tax_rank), aldex2_res, envir = .GlobalEnv)
  
  # Ensure required columns exist
  if (!all(c("ef_aldex", "padj", "feature") %in% colnames(aldex2_res))) {
    message("Missing expected columns in ALDEx2 result.")
    return(invisible(NULL))
  }
  
  # Prepare results for plotting
  df_fig <- aldex2_res %>%
    mutate(
      ef_aldex = as.numeric(ef_aldex),
      padj = as.numeric(padj),
      Enrichment = ifelse(ef_aldex > 0, "Plaque", "Abscess")
    )
  
  if (nrow(df_fig) == 0) {
    message("No significant taxa found in ALDEx2 results.")
    return(invisible(NULL))
  }
  
  df_fig$feature <- factor(df_fig$feature, levels = df_fig$feature)
  df_fig$Enrichment <- factor(df_fig$Enrichment, levels = c("Plaque", "Abscess"))
  
  # Dot plot
  dot_plot <- df_fig %>%
    ggplot(aes(x = abs(ef_aldex), y = feature)) +
    geom_point(aes(size = -log10(padj), color = Enrichment)) +
    scale_color_manual(values = plot_colors) +
    scale_size(range = c(1, 5), name = "-log10(P-value)") +
    labs(title = NULL, x = "EF Aldex", y = NULL, color = "Enrichment") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Bar plot
  bar_plot <- df_fig %>%
    ggplot(aes(x = feature, y = abs(ef_aldex), fill = Enrichment)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +
    coord_flip() +
    labs(x = NULL, y = "EF Aldex", title = NULL) +
    scale_fill_manual(values = plot_colors) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  #Return
  return(list(results = df_fig, bar_plot = bar_plot, dot_plot = dot_plot))
}


## ---- Combine DA Methods----

combine_DA <- function(maaslin2_results, ancombc2_results, aldex2_results, group, tax_level, filter_confidence = NULL){
  
  #Function for cleaning the taxa in the various result tables to keep things consistent 
  clean_taxon <- function(x) {
    x  %>%
      stri_trans_nfkc() %>%                 # Unicode normalize
      stringr::str_replace_all("^(g_|s_|p_)", "") %>%
      str_replace_all("\\p{Z}+", " ")  %>%  # any Unicode space class -> single space
      str_squish() %>%                     # trim + collapse internal spaces
      str_replace_all("[\\p{P}\\p{S}]+", " ")  %>%  # drop punctuation/symbols
      str_replace_all("-", "")  %>%        
      str_replace_all(" s ", " ")  %>%
      str_replace_all(" F ", " F")  %>%
      str_replace_all(" G ", " G")  %>%
      str_squish() 
  }
  
    # Safe limit function
    safe_lim <- function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return(0)
      max(abs(x))
    }
    
    # Helper: p-value to stars
    p_to_stars <- function(p) {
      ifelse(is.na(p), "",
             ifelse(p <= 0.001, "***",
                    ifelse(p <= 0.01,  "**",
                           ifelse(p <= 0.05, "*", ""))))
    }

    # Prepare the two dataframes for merging
    # MAASLIN2
    maaslin2_results_clean <- maaslin2_results %>%
      dplyr::rename(taxon = feature) %>%
      dplyr::mutate(taxon = clean_taxon(taxon)) %>%
      dplyr::select(taxon,
             Maaslin2_coef = coef,
             Maaslin2_pval = pval,
             Maaslin2_qval = qval,
             Maaslin2_value = value)
    
    # ANCOMBC2  (use all_of for the dynamic names!)
    dff_col_name <- grep(paste0("^diff_", group), colnames(ancombc2_results), value = TRUE)
    lfc_col_name <- grep(paste0("^lfc_",  group), colnames(ancombc2_results), value = TRUE)
    ss_col_name  <- grep(paste0("^passed_ss_", group), colnames(ancombc2_results), value = TRUE)
    p_col_name   <- grep(paste0("^p_",    group), colnames(ancombc2_results), value = TRUE)
    q_col_name   <- grep(paste0("^q_",    group), colnames(ancombc2_results), value = TRUE)
    
    ancombc2_results_clean <- ancombc2_results %>%
      dplyr::mutate(taxon = clean_taxon(taxon)) %>%
      dplyr::select(taxon,
             ANCOMBC2_DFF       = all_of(dff_col_name),
             ANCOMBC2_LFC       = all_of(lfc_col_name),
             ANCOMBC2_p_val     = all_of(p_col_name),
             ANCOMBC2_q_val     = all_of(q_col_name),
             ANCOMBC2_passed_SS = all_of(ss_col_name)) %>%
      dplyr::filter(ANCOMBC2_DFF %in% TRUE) 
    
    # ALDEx2
    aldex2_res_clean <- aldex2_results %>%
      dplyr::mutate(ef_aldex = as.numeric(ef_aldex),
             padj     = as.numeric(padj),
             feature  = as.character(feature)) %>%
      dplyr::rename(taxon = feature) %>%
      dplyr::mutate(taxon = clean_taxon(taxon),
             Enrichment = if_else(ef_aldex > 0, "Plaque", "Abscess")) %>%
      dplyr::select(taxon, Aldex2_Enrichment = Enrichment, Aldex2_EF = ef_aldex, Aldex2_padj = padj)
      
    
    # ---- Join all of them together!
    
      first_non_na <- function(x) { x[which(!is.na(x))[1]] %||% NA }
      
      DA_results_df <- dplyr::full_join(maaslin2_results_clean, ancombc2_results_clean, by = "taxon") %>%
        dplyr::full_join(aldex2_res_clean, by = "taxon") %>%
        group_by(taxon) %>%
        dplyr::summarise(across(everything(), first_non_na), .groups = "drop") %>%
        dplyr::mutate(confidence = case_when(
          #We are counting NAs, so if there are 2 NAs it's low confidence, and 0 it's high
          rowSums(is.na(pick(Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) == 2 ~ "Low",
          rowSums(is.na(pick(Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) == 1 ~ "Medium",
          rowSums(is.na(pick(Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) == 0 ~ "High"
        ))

    DA_results_df <- DA_results_df %>% dplyr::rowwise() %>%
      dplyr::mutate(
        all_positive_or_negative = {
          vals <- dplyr::c_across(c(Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))
          vals <- vals[!is.na(vals)]
          length(vals) > 0 && (all(vals > 0) || all(vals < 0))
        }
      ) %>%
      ungroup()
    
    # Warning just in case the directions of the tests don't agree
    if (any(DA_results_df$all_positive_or_negative == FALSE, na.rm = TRUE)) {
      warning("Some values in 'all_positive_or_negative' are FALSE!")}
    
    # Assign the results to the environment
    assign(paste0("DA_results_", tax_level), DA_results_df, envir = .GlobalEnv)
    
    
    ## --- Plot the confidence scores of all of them

    # Only filter if filter_confidence is provided
    df <- DA_results_df
    if (!is.null(filter_confidence)) {
      df <- df %>% filter(!confidence %in% filter_confidence)
    }

    
    # Filter if requested
    df <- DA_results_df
    if (!is.null(filter_confidence)) {
      df <- df %>% filter(!confidence %in% filter_confidence)
    }
    
    # Prepare data
    df_plot <- df %>%
      dplyr::mutate(
        stars_maaslin  = p_to_stars(Maaslin2_qval),
        stars_ancombc2 = p_to_stars(ANCOMBC2_q_val),
        stars_aldex2   = p_to_stars(Aldex2_padj),
        enrichment = case_when(
          Maaslin2_coef > 0 | Aldex2_EF > 0 | ANCOMBC2_LFC > 0 ~ "Plaque",
          Maaslin2_coef < 0 | Aldex2_EF < 0 | ANCOMBC2_LFC < 0 ~ "Abscess",
          TRUE ~ NA_character_
        ),
        confidence = factor(confidence, levels = c("High","Medium","Low"))
      )
    
    # Column positions / labels
    x_pos <- c(Maaslin2 = 1, ANCOMBC2 = 2, ALDEx2 = 3)
    x_lab <- c("MaAsLin2", "ANCOMBC2", "ALDEx2")
    
    
    
    ## --- Set pal limits before plotting --- 
    
    # Subset ab only
    df_plot_abscess <- df_plot %>% filter(enrichment == "Abscess")
    
    
    #Switch the sign (since all the abscess values are negative)
    df_plot_abscess <- df_plot_abscess %>% mutate(
      Maaslin2_coef = Maaslin2_coef * -1 ,
      ANCOMBC2_LFC = ANCOMBC2_LFC * -1 ,
      Aldex2_EF = Aldex2_EF * -1,
    )
    
    
    # Subset plaque only
    df_plot_plaque <- df_plot %>% filter(enrichment == "Plaque")
    
    # Combine plaque + abscess to get shared limits
    lim_maas <- safe_lim(c(df_plot_plaque$Maaslin2_coef, df_plot_abscess$Maaslin2_coef))
    lim_ancb <- safe_lim(c(df_plot_plaque$ANCOMBC2_LFC,  df_plot_abscess$ANCOMBC2_LFC))
    lim_aldx <- safe_lim(c(df_plot_plaque$Aldex2_EF,     df_plot_abscess$Aldex2_EF))
    
    ## --- Plaque Plot --- ##
    
    pal <- list(
      low = "#ddebfe", 
      high = "#3185FC",
      na = "white"
    )

    # Exact global order used by bubbles (per-enrichment), reversed once
    y_levels_plaque <- df_plot %>%
      dplyr::filter(enrichment == "Plaque") %>%
      dplyr::arrange(confidence, taxon) %>%   # alphabetical within each confidence
      dplyr::pull(taxon) %>%
      unique()
    

    df_plot_plaque <- df_plot_plaque %>%
      dplyr::mutate(taxon = forcats::fct_rev(factor(taxon, levels = y_levels_plaque)))

    print(df_plot_plaque %>% filter(confidence == "Low"))
    
    
     #I guess you have to do this to get the ordering correct! So annoying.
    scaffold_plaque <- df_plot_plaque %>% distinct(confidence, taxon)
    
    
    # Build plot
    p_plaque <- ggplot() +
      geom_blank(data = scaffold_plaque, aes(x = 0, y = taxon)) +  # again, for ordering

      # Column 1: MaAsLin2
      geom_tile(
        data = df_plot_plaque %>% filter(!is.na(Maaslin2_coef)),
        aes(x = x_pos["Maaslin2"], y = taxon, fill = Maaslin2_coef),
        color = "black",    # border color
        width = 1,            # <— make the tile fill the column
        height = 1, 
        size = 0.3              # border thickness
      ) +
      scale_fill_gradient2(
        name = "MaAsLin2",
        low = pal$low,
        high = pal$high,
        na.value = pal$na,
        limits = c(0, lim_maas)
      ) +
      geom_text(
        data = df_plot_plaque,
        aes(x = x_pos["Maaslin2"], y = taxon, label = stars_maaslin),
        size = 5, vjust = .8
      ) +
      
      # Column 2: ANCOMBC2
      new_scale_fill() +
      geom_tile(
        data = df_plot_plaque %>% filter(!is.na(ANCOMBC2_LFC)),
        aes(x = x_pos["ANCOMBC2"], y = taxon, fill = ANCOMBC2_LFC),
        color = "black",    # border color
        width = 1,            # <— make the tile fill the column
        height = 1, 
        size = 0.3              # border thickness
      ) +
      scale_fill_gradient2(
        name = "ANCOMBC2",
        low = pal$low,
        high = pal$high,
        na.value = pal$na,
        limits = c(0, lim_ancb)
      ) +
      geom_text(
        data = df_plot_plaque,
        aes(x = x_pos["ANCOMBC2"], y = taxon, label = stars_ancombc2),
        size = 5, vjust = .8
      ) +
      
      # Column 3: ALDEx2
      new_scale_fill() +
      geom_tile(
        data = df_plot_plaque %>% filter(!is.na(Aldex2_EF)),
        aes(x = x_pos["ALDEx2"], y = taxon, fill = Aldex2_EF),
        color = "black",    # border color
        width = 1,            # <— make the tile fill the column
        height = 1, 
        size = 0.3              # border thickness
      ) +
      scale_fill_gradient2(
        name = "ALDEx2",
        low = pal$low,
        high = pal$high,
        na.value = pal$na,
        limits = c(0, lim_aldx)
      ) +
      geom_text(
        data = df_plot_plaque,
        aes(x = x_pos["ALDEx2"], y = taxon, label = stars_aldex2),
        size = 5, vjust = .8
      ) +
      scale_x_continuous(
        breaks = x_pos, labels = x_lab,
        limits = c(min(x_pos) - 0.5, max(x_pos) + 0.5),  # <— half-step limits
        expand = expansion(mult = c(0, 0))               # <— no padding
      ) +
      scale_y_discrete(
        drop   = TRUE,
        expand = expansion(mult = c(0, 0))  # <-- add this
      ) +
      facet_grid(confidence ~ ., scales = "free_y", space = "free_y") +
      labs(x = NULL, y = NULL) +
      coord_cartesian(clip = "off") +
      theme_bw(base_size = 12) +
      theme(
        text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),   # removes the rectangle border
        panel.spacing.x = unit(0.6, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(face = "italic", color = "black", size = 12)) + 
      theme(
        legend.position = "bottom",
        legend.box = "vertical"  
      )
    

    ## --- Abscess Plot --- ##
  
    pal <- list(
      low = "#ffe6e9", 
      high = "#FF495C",
      na = "white"
    )
  

    # Exact global order used by bubbles (per-enrichment), reversed once
    y_levels_abscess <- df_plot %>%
      dplyr::filter(enrichment == "Abscess") %>%
      dplyr::arrange(confidence, taxon) %>%   # alphabetical within each confidence
      dplyr::pull(taxon) %>%
      unique()
    
    df_plot_abscess <- df_plot_abscess %>%
      dplyr::mutate(taxon = forcats::fct_rev(factor(taxon, levels = y_levels_abscess)))
    
    #I guess you have to do this to get the ordering correct! So annoying.
    scaffold_abscess <- df_plot_abscess %>% distinct(confidence, taxon)
    
    # Build plot
    a_plot <- ggplot() +
      geom_blank(data = scaffold_abscess, aes(x = 0, y = taxon)) +  # again, for ordering
      
      # Column 1: MaAsLin2
      geom_tile(
        data = df_plot_abscess %>% filter(!is.na(Maaslin2_coef)),
        aes(x = x_pos["Maaslin2"], y = taxon, fill = Maaslin2_coef),
        color = "black",    # border color
        width = 1,            # <— make the tile fill the column
        height = 1, 
        size = 0.3              # border thickness
      ) +
      scale_fill_gradient2(
        name = "MaAsLin2",
        low = pal$low,
        high = pal$high,
        na.value = pal$na,
        limits = c(0, lim_maas)
      ) +
      geom_text(
        data = df_plot_abscess,
        aes(x = x_pos["Maaslin2"], y = taxon, label = stars_maaslin),
        size = 5, vjust = .8
      ) +
      
      # Column 2: ANCOMBC2
      new_scale_fill() +
      geom_tile(
        data = df_plot_abscess %>% filter(!is.na(ANCOMBC2_LFC)),
        aes(x = x_pos["ANCOMBC2"], y = taxon, fill = ANCOMBC2_LFC),
        color = "black",    # border color
        width = 1,            # <— make the tile fill the column
        height = 1, 
        size = 0.3              # border thickness
      ) +
      scale_fill_gradient2(
        name = "ANCOMBC2",
        low = pal$low,
        high = pal$high,
        na.value = pal$na,
        limits = c(0, lim_ancb)
      ) +
      geom_text(
        data = df_plot_abscess,
        aes(x = x_pos["ANCOMBC2"], y = taxon, label = stars_ancombc2),
        size = 5, vjust = .8
      ) +
      
      # Column 3: ALDEx2
      new_scale_fill() +
      geom_tile(
        data = df_plot_abscess %>% filter(!is.na(Aldex2_EF)),
        aes(x = x_pos["ALDEx2"], y = taxon, fill = Aldex2_EF),
        color = "black",    # border color
        width = 1,            # <— make the tile fill the column
        height = 1, 
        size = 0.3              # border thickness
      ) +
      scale_fill_gradient2(
        name = "ALDEx2",
        low = pal$low,
        high = pal$high,
        na.value = pal$na,
        limits = c(0, lim_aldx)
      ) +
      geom_text(
        data = df_plot_abscess,
        aes(x = x_pos["ALDEx2"], y = taxon, label = stars_aldex2),
        size = 5, vjust = .8
      ) +
      scale_x_continuous(
        breaks = x_pos, labels = x_lab,
        limits = c(min(x_pos) - 0.5, max(x_pos) + 0.5),  # <— half-step limits
        expand = expansion(mult = c(0, 0))               # <— no padding
      ) +
      scale_y_discrete(
        drop   = TRUE,
        expand = expansion(mult = c(0, 0))  # <-- add this
      ) +
  facet_grid(confidence ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),   # removes the rectangle border
    panel.spacing.x = unit(0.6, "lines"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(face = "italic", color = "black", size = 12)
  ) + 
  theme(
    legend.position = "bottom",
    legend.box = "vertical"  
  )

    
    return(list(results = DA_results_df, plaque_plot = p_plaque, abscess_plot = a_plot, df_plot = df_plot))
}

# This is the function I actually used for plotting - combines the results of the function above with a bubble plot for asessing robustness.
combine_DA_with_bubbles <- function(da, 
    CSV_formatted_relab_df,
    group, tax_level,
    filter_confidence = NULL,
    heatmap_side = c("Plaque","Abscess"),
    bubble_max_size = 8,
    plot_layout_widths = c(2, 1)){
  
  
  require(dplyr); require(stringr); require(tidyr)
  require(ggplot2); require(forcats); require(patchwork)
  require(ggnewscale); require(purrr)
  
  heatmap_side <- match.arg(heatmap_side)
  
  # ------ Get output of combine_DA function
  heatmap_plot <- if (heatmap_side == "Plaque") da$plaque_plot else da$abscess_plot
  
  df_plot <- da$df_plot %>%
    mutate(
      confidence = factor(confidence, levels = c("High","Medium","Low")),
      taxon = as.character(taxon)
    )
  
  # --------- Bubble plot data with robust CANONICAL MAPPING (fixes weird names)
  # One canonical normalizer used on BOTH sides
  canon_taxon <- function(x) {
    x %>%
      stringi::stri_trans_nfkc() %>%          # Unicode normalize
      stringr::str_to_lower() %>%             # casefold
      # strip rank prefixes (at start or after ; | or space)
      stringr::str_replace_all("(^|[;|\\|[:space:]])(?:[kpcofgs]__)", "\\1") %>%
      # legacy single-underscore prefixes
      stringr::str_replace_all("^(g_|s_|p_)", "") %>%
      # species suffix artifact
      stringr::str_replace_all("_s__", " ") %>%
      # TM7 / G group variants: normalize separators like G1, G-1, G 1 -> g1
      stringr::str_replace_all("\\b[gG][\\s_-]*([0-9]+)\\b", " g\\1") %>%
      # normalize tm7 tokens (tm-7, tm_7, tm 7) -> tm7
      stringr::str_replace_all("\\btm[\\s_-]*([0-9]+)\\b", " tm\\1") %>%
      # remove brackets
      stringr::str_replace_all("\\[|\\]", "") %>%
      # underscores/hyphens -> space; drop symbols/punct
      stringr::str_replace_all("[_\\-]+", " ") %>%
      stringr::str_replace_all("[\\p{P}\\p{S}]+", " ") %>%
      stringr::str_squish()
  }
  
  enrich_pick <- if (heatmap_side == "Plaque") "Plaque" else "Abscess"
  
  # DA taxon set (only those shown on the selected side), with canonical key
  da_taxa <- df_plot %>%
    filter(enrichment == enrich_pick) %>%
    distinct(taxon, confidence) %>%
    mutate(taxon_canon = canon_taxon(taxon))
  
  # --- Build a column map for the CSV taxa *first*, so we can collapse synonyms
  tax_cols <- grep("^(k|p|c|o|f|g|s)__|^(g_|s_|p_)", names(CSV_formatted_relab_df), value = TRUE)
  
  csv_tax_map <- tibble(
    tax_orig  = tax_cols,
    tax_canon = canon_taxon(tax_cols)
  )
  
  # Keep only CSV taxa whose canonical key appears in DA on this side
  keep_keys <- unique(da_taxa$taxon_canon)
  csv_tax_map_keep <- csv_tax_map %>% filter(tax_canon %in% keep_keys)
  
  
  # --- Melt CSV, canonicalize taxon names, standardize Type, collapse synonyms
  csv_long <- CSV_formatted_relab_df %>%
    pivot_longer(cols = all_of(csv_tax_map_keep$tax_orig),
                 names_to = "tax_orig", values_to = "relab") %>%
    dplyr::left_join(csv_tax_map_keep, by = "tax_orig") %>%
    mutate(
      # Standardize Type values robustly
      Type = case_when(
        tolower(Type) %in% c("p","plaque")  ~ "Plaque",
        tolower(Type) %in% c("a","abscess") ~ "Abscess",
        TRUE ~ as.character(Type)
      ),
      Sample = str_replace_all(Sample, "-P|-[Aa]", "")
    ) %>%
    group_by(Sample, Type, tax_canon) %>%
    dplyr::summarise(relab = sum(relab, na.rm = TRUE), .groups = "drop")  # collapse synonyms -> one key
  
  # Pivot to wide by Type
  csv_wide <- csv_long %>%
    pivot_wider(names_from = "Type", values_from = "relab", values_fill = 0)
  
  # Ensure both Abscess and Plaque columns exist
  if (!"Abscess" %in% names(csv_wide)) csv_wide$Abscess <- 0
  if (!"Plaque"  %in% names(csv_wide)) csv_wide$Plaque  <- 0
  
  data_pull <- csv_wide %>%
    mutate(
      Ratio     = ifelse(Plaque == 0 & Abscess == 0, NA_real_, Abscess / dplyr::na_if(Plaque, 0)),
      RatioFlag = case_when(
        is.na(Ratio) ~ NA_integer_,
        Ratio > 1    ~ 1L,   # Abscess enriched
        TRUE         ~ -1L    # Plaque enriched (ties & <=1)
      )
    )
  
  df_summary_canon <- data_pull %>%
    dplyr::group_by(tax_canon) %>%
    dplyr::summarise(
      Abscess_enriched = sum(RatioFlag == 1,  na.rm = TRUE),
      Plaque_enriched  = sum(RatioFlag == -1, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join back to DA (original taxon label + confidence) by canonical key
  bubble_df <- da_taxa %>%
    dplyr::left_join(df_summary_canon, by = c("taxon_canon" = "tax_canon")) %>%
    dplyr::mutate(
      Abscess_enriched = replace_na(Abscess_enriched, 0L),
      Plaque_enriched  = replace_na(Plaque_enriched, 0L)
    ) %>%
    select(taxon, confidence, Abscess_enriched, Plaque_enriched) %>%
    pivot_longer(
      cols = c(Abscess_enriched, Plaque_enriched),
      names_to = "group2", values_to = "count"
    ) %>%
    mutate(
      group2     = recode(group2, "Abscess_enriched"="A", "Plaque_enriched"="P"),
      confidence = factor(confidence, levels = c("High","Medium","Low"))
    )
  

  # -------- Exact y-order, and per-facet scaffold (no extras; height matches)
  y_levels <- df_plot %>%
    filter(enrichment == enrich_pick) %>%
    arrange(enrichment, confidence, taxon) %>%
    pull(taxon) %>% unique() %>% rev()
  
  # Lock bubble_df to global levels (order identical to heatmap)
  bubble_df <- bubble_df %>%
    mutate(taxon = factor(taxon, levels = y_levels))
  
  taxa_by_conf <- df_plot %>%
    dplyr::filter(enrichment == enrich_pick) %>%
    dplyr::group_by(confidence) %>%
    dplyr::summarise(
      facet_taxa = list(y_levels[y_levels %in% unique(taxon)]),
      .groups = "drop"
    )
  
  scaffold <- map_dfr(seq_len(nrow(taxa_by_conf)), function(i) {
    conf_i <- taxa_by_conf$confidence[i]
    tax_i  <- taxa_by_conf$facet_taxa[[i]]
    expand.grid(
      confidence = factor(conf_i, levels = c("High","Medium","Low")),
      taxon      = factor(tax_i, levels = y_levels),
      group2     = c("A","P"),
      stringsAsFactors = FALSE
    )
  })
  
  bubble_df_full <- scaffold %>%
    dplyr::left_join(bubble_df, by = c("confidence","taxon","group2")) %>%
    dplyr::mutate(count = replace_na(count, 0L))
  
  bubble_A <- bubble_df_full %>% filter(group2 == "A")
  bubble_P <- bubble_df_full %>% filter(group2 == "P")
  

  # ------- Bubble plot (independent color ramps; shared size)
  bubble_plot <- ggplot() +
    # Abscess (red)
    geom_point(
      data = bubble_A,
      aes(x = group2, y = taxon, size = count, fill = count),
      shape = 21, alpha = 0.9, stroke = 0.4, color = "black"
    ) +
    # geom_text(
    #   data = bubble_A,
    #   aes(x = group2, y = taxon, label = count),
    #   size = 3, vjust = 0.5
    # ) +
    scale_fill_gradientn(
      colors = c("white", "#FF495C"),
      limits = c(0, 25),
      breaks = seq(0, 25, 5),        # <-- force 25 to appear
      oob = scales::squish,
      guide = "none"
    ) +
    ggnewscale::new_scale("fill") +
    # Plaque (blue)
    geom_point(
      data = bubble_P,
      aes(x = group2, y = taxon, size = count, fill = count),
      shape = 21, alpha = 0.9, stroke = 0.4, color = "black"
    ) +
    # geom_text(
    #   data = bubble_P,
    #   aes(x = group2, y = taxon, label = count),
    #   size = 3, vjust = 0.5
    # ) +
    scale_fill_gradientn(
      colors = c("white", "#3185FC"),
      limits = c(0, 25),
      breaks = seq(0, 25, 5),        # <-- force 25 to appear
      oob = scales::squish,
      guide = "none"
    ) +
    scale_size_area(max_size = bubble_max_size, limits = c(0, NA)) +
    scale_x_discrete(labels = c(A = "A", P = "P")) +
    facet_grid(confidence ~ ., scales = "free_y", space = "free_y") +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, size = "Count") +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(face = "italic", size = 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      panel.border = element_blank()
    )
  
  
  # Size legend (make sure it shows 25)
  scale_size_area(
    name = "Count",
    max_size = bubble_max_size,
    limits = c(0, 25),
    breaks = seq(0, 25, 5),        # <-- include 25 explicitly
    guide = guide_legend(
      override.aes = list(fill = "grey95", color = "black")
    )
  )
  
  bubble_plot <- bubble_plot +
    scale_size_area(name="Count", limits=c(0,25), breaks=c(0,5,10,15,25),
                    guide = guide_legend(nrow = 3, byrow = TRUE, direction = "horizontal"))
  
  #Have to mask the heatmap legend so I can actually see the bubble plot legend and use it for scaling
  heatmap_masklegend <- heatmap_plot +
    theme(legend.position = "none")
  
  combined <- (heatmap_masklegend + bubble_plot + plot_layout(widths = plot_layout_widths, guides = "collect")) &
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  heatmap_noleg <- heatmap_plot + theme(legend.position = "none", 
                                        axis.title.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.text.y  = element_blank(),
                                        axis.text.x  = element_blank())
  bubble_noleg  <- bubble_plot  + theme(legend.position = "none", 
                                        axis.title.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.text.y  = element_blank())
  
  combined_nolabs <- heatmap_noleg + 
    bubble_noleg 

  list(
    results       = da$results,
    heatmap       = heatmap_plot,
    bubble_plot   = bubble_plot,
    combined_plot = combined,
    combined_plot_nolabs = combined_nolabs,
    bubble_data   = bubble_df_full
  )
}




# ---- Random Forest ----
rank_rf_feature_table <- function(df, folds, mtry = NULL, top_n = NULL, add_mean = TRUE) {
  
  if (!is.factor(df$Type)) df$Type <- factor(df$Type)
  p <- ncol(df) - 1L
  if (is.null(mtry)) mtry <- max(1L, floor(sqrt(p)))
  
  
  get_fold_ranks <- function(train_idx, fold_id) {
    fit <- caret::train(
      Type ~ .,
      data = df[train_idx, ],
      method = "rf",
      trControl = caret::trainControl(method = "none"),
      tuneGrid = data.frame(mtry = mtry),
      importance = TRUE
    )
    imp <- caret::varImp(fit)$importance %>%
      rownames_to_column("Feature")
    
    # collapse to a single score column named Score (usually "Overall" already)
    score_cols <- setdiff(names(imp), "Feature")
    if (length(score_cols) == 0) stop("No importance scores returned.")
    if (length(score_cols) == 1) {
      imp <- transmute(imp, Feature, Score = .data[[score_cols]])
    } else {
      imp <- imp %>% mutate(Score = rowMeans(across(all_of(score_cols)), na.rm = TRUE)) %>%
        select(Feature, Score)
    }
    
    # rank (1 = most important)
    imp <- imp[order(imp$Score, decreasing = TRUE), , drop = FALSE]
    imp$Rank <- seq_len(nrow(imp))
    
    # return Feature + FoldX rank column
    fold_col <- paste0("Fold", fold_id)
    out <- imp[, c("Feature", "Rank")]
    names(out)[2] <- fold_col
    out
  }
  
  fold_tables <- Map(function(idx, i) get_fold_ranks(idx, i), folds, seq_along(folds))
  
  # Join by Feature (NOT by rank) so each taxon is one row
  wide <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE), fold_tables)
  
  # Optionally keep only taxa that are in top_n in at least one fold
  if (!is.null(top_n)) {
    rank_cols <- grep("^Fold\\d+$", names(wide), value = TRUE)
    keep <- apply(wide[rank_cols], 1, function(v) any(!is.na(v) & v <= top_n))
    wide <- wide[keep, , drop = FALSE]
  }
  
  # Mean rank across folds (lower is better)
  if (add_mean) {
    rank_cols <- grep("^Fold\\d+$", names(wide), value = TRUE)
    wide$MeanRank <- rowMeans(wide[rank_cols], na.rm = TRUE)
    wide <- wide[order(wide$MeanRank), , drop = FALSE]
  }
  
  rownames(wide) <- NULL
  wide
}


# ---- Topic Modeling ----

#This function helps you optimize the number of topics (k) and also the scaling factor by which you multiply your input your data (essentially the threshold). 
iterate_scaling_factors <- function(phy_obj, 
                                    taxa_level, 
                                    scale_factors, #A list of data scaling factors. The higher the values the longer this will take.
                                    k_values, #A list of topic numbers to use. You need this regardless of setting auto_chose_k to TRUE or FALSE
                                    auto_chose_k = TRUE, #Set this to TRUE if you want to automatically chose the number of topics (k). If you want to iterate over k_values and then have the model computed for each one, then chose FALSE.
                                    type_column, 
                                    method = "Gibbs") {
  results_list <- list()
  
  library(cluster)
  library(Matrix)
  library(R.utils)
  
  # Extract metadata
  meta_data <- phy_obj@sam_data %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(c(type_column, "Sample")))
  
  base_counts_data <- tax_glom(phy_obj, taxa_level) %>%
    norm_tss() %>%
    psmelt() %>%
    dplyr::select(c("Sample", "Abundance", all_of(taxa_level))) %>%
    pivot_wider(names_from = all_of(taxa_level), values_from = "Abundance") %>%
    column_to_rownames(var = "Sample")
  
  for (scaling_factor in scale_factors) {
    
    # Create Relative Abundance Table
    counts_data <- base_counts_data%>%
      mutate(across(where(is.numeric), ~ round(. * scaling_factor)))
    
    # Store and reset rownames
    row_names <- rownames(counts_data)
    counts_data <- as.data.frame(lapply(counts_data, as.integer))
    rownames(counts_data) <- row_names
    
    
    # ------------ If we want to auto chose k, we are going to use the following code to automatically chose the best k-value and then feed that into the model
    if (isTRUE(auto_chose_k))  {
      
      cat("Auo-chosing topics number using a scaling factor of", scaling_factor, "\n")
      
      withTimeout({
        #Throw an error if the model takes way too long.
        result <- FindTopicsNumber(
          counts_data,
          topics = k_values,
          metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
          method = "Gibbs",
          control = list(seed = 1234),
          mc.cores = 2L,
          verbose = TRUE
        )}, timeout = 120, onTimeout = "error") 
      
      # Combine each metric into a named row of a data frame
      topic_numbers <- result$topics
      
      df_metrics <- data.frame(
        Metric = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
        stringsAsFactors = FALSE
      )
      
      # Extract each metric as a named row
      df_metrics_values <- rbind(
        setNames(as.data.frame(t(result$Griffiths2004)), paste0("k_", topic_numbers)),
        setNames(as.data.frame(t(result$CaoJuan2009)), paste0("k_",topic_numbers)),
        setNames(as.data.frame(t(result$Arun2010)), paste0("k_",topic_numbers)),
        setNames(as.data.frame(t(result$Deveaud2014)), paste0("k_",topic_numbers))
      )
      
      # Combine metric names and values
      df_result <- cbind(df_metrics, df_metrics_values)
      
      #Find optimial topic number where:
      ##Griffiths2004 and Deveaud2014 are maximized
      ##CaoJuan2009 and Arun2010 are minimized
      
      #Transpose df_result to long format
      df_long <- df_result %>%
        pivot_longer(cols = starts_with("k"), names_to = "k", values_to = "score")
      
      df_long %>%
        group_by(Metric) 
      
      # Normalize each metric (min-max scaling)
      df_scaled <- df_long %>%
        group_by(Metric) %>%
        mutate(score_scaled = (score - min(score)) / (max(score) - min(score))) %>%
        ungroup()
      
      # Flip scores for metrics to minimize
      df_scaled <- df_scaled %>%
        mutate(
          score_final = case_when(
            Metric %in% c("CaoJuan2009", "Arun2010") ~ 1 - score_scaled,
            TRUE ~ score_scaled
          )
        )
      
      # Average score per k across all metrics
      df_combined <- df_scaled %>%
        group_by(k) %>%
        dplyr::summarise(
          combined_score = mean(score_final),
          .groups = "drop"
        ) %>%
        arrange(dplyr::desc(combined_score))
      
      # Best topic number overall
      best_k <- df_combined$k[1]
      k_value <- as.numeric(gsub("k_", "", best_k))
      
      cat("And the best topic number is....", k_value, "!", "\n")
      
      # ------- Now that we've found the best k-number, we can do all the other fun stuff by plugging that into the model
      # Create topic model
      counts_mat <- as.matrix(counts_data)  # Convert to regular dense matrix
      dtm <- as(counts_mat, "dgCMatrix")    # Convert to sparse matrix
      withTimeout({
        #Throw an error if the model takes way too long.
        lda_model <- topicmodels::LDA(dtm, k = k_value, method = method, control = list(seed = 1234))
      }, timeout = 120, onTimeout = "error")
      
      perplex <- topicmodels::perplexity(lda_model, newdata = as.matrix(counts_data))
      
      # Extract gamma matrix (topic proportions)
      g_df <- data.frame(tidytext::tidy(lda_model, matrix = "gamma")) %>%
        arrange(document, topic) %>%
        pivot_wider(names_from = topic, values_from = gamma)
      
      # Merge topic proportions with metadata
      res_with_annotations <- meta_data %>%
        rownames_to_column(var="document") %>%
        merge(g_df, by = "document")
      
      topic_matrix <- res_with_annotations %>%
        dplyr::select(-all_of(type_column)) %>%
        dplyr::select(-Sample) %>%
        column_to_rownames("document")
      
      # Silhouette score calculation
      dist_mat <- dist(topic_matrix)
      cluster_labels <- factor(res_with_annotations[[type_column]])
      sil <- silhouette(as.numeric(cluster_labels), dist_mat)
      sil_df <- as.data.frame(sil)
      sil_df$Type <- cluster_labels[as.numeric(rownames(sil_df))]
      
      # Median silhouette per Type
      # Rename the silhouette width column properly
      colnames(sil_df)[colnames(sil_df) == "sil_width" | colnames(sil_df) == "silhouette.width"] <- "silhouette"
      
      sil_medians <- sil_df %>%
        dplyr::group_by(Type) %>%
        dplyr::summarise(
          median_silhouette = median(silhouette, na.rm = TRUE),
          .groups = "drop") %>% 
        tidyr::pivot_wider(names_from = Type, values_from = median_silhouette, names_prefix = "sil_median_") %>%
        mutate(
          k_value = k_value,
          scaling_factor = scaling_factor,
          perplexity = perplex
        )
      
      cat("Model created, moving to next. \n")
      
      # Store results
      results_list[[paste0("k", k_value, "_scale", scaling_factor)]] <- sil_medians
      
      # ------------ If you don't want to automatically chose K, this forces the model to build with each number in k-values
    } else {
      for (k_value in k_values) {
        
        cat("running model with ", k_value, " topics and scaling factor of ", scaling_factor, "\n")
        
        # Create topic model
        counts_mat <- as.matrix(counts_data)  # Convert to regular dense matrix
        dtm <- as(counts_mat, "dgCMatrix")    # Convert to sparse matrix
        withTimeout({
          #Throw an error if the model takes way too long.
          lda_model <- topicmodels::LDA(dtm, k = k_value, method = method, control = list(seed = 1234))
        }, timeout = 120, onTimeout = "error")
        
        perplex <- topicmodels::perplexity(lda_model, newdata = as.matrix(counts_data))
        
        # Extract gamma matrix (topic proportions)
        g_df <- data.frame(tidytext::tidy(lda_model, matrix = "gamma")) %>%
          arrange(document, topic) %>%
          pivot_wider(names_from = topic, values_from = gamma)
        
        # Merge topic proportions with metadata
        res_with_annotations <- meta_data %>%
          rownames_to_column(var="document") %>%
          merge(g_df, by = "document")
        
        topic_matrix <- res_with_annotations %>%
          dplyr::select(-all_of(type_column)) %>%
          dplyr::select(-Sample) %>%
          column_to_rownames("document")
        
        # Silhouette score calculation
        dist_mat <- dist(topic_matrix)
        cluster_labels <- factor(res_with_annotations[[type_column]])
        sil <- silhouette(as.numeric(cluster_labels), dist_mat)
        sil_df <- as.data.frame(sil)
        sil_df$Type <- cluster_labels[as.numeric(rownames(sil_df))]
        
        # Median silhouette per Type
        # Rename the silhouette width column properly
        colnames(sil_df)[colnames(sil_df) == "sil_width" | colnames(sil_df) == "silhouette.width"] <- "silhouette"
        
        sil_medians <- sil_df %>%
          dplyr::group_by(Type) %>%
          dplyr::summarise(
            median_silhouette = median(silhouette, na.rm = TRUE),
            .groups = "drop") %>% 
          tidyr::pivot_wider(names_from = Type, values_from = median_silhouette, names_prefix = "sil_median_") %>%
          mutate(
            k_value = k_value,
            scaling_factor = scaling_factor,
            perplexity = perplex
          )
        
        # Store results
        results_list[[paste0("k", k_value, "_scale", scaling_factor)]] <- sil_medians
      }
    } 
  }
  
  # Combine all results
  final_df <- bind_rows(results_list)
  return(list(results = final_df))
}


##Prep data and convert to whole number by factor
#Prepare data from phyloseq object and binarize at a certain threshold of relative abundance
prep_data_scale <- function(phy_obj, taxa_level, scaling_factor, type_column){
  # Extract metadata
  meta_data <- phy_obj@sam_data %>% 
    as.matrix() %>% as.data.frame() %>% 
    dplyr::select(c(type_column, "Sample"))

  #Create Relative Abundance Table
  counts_data <- tax_glom(phy_obj, taxa_level) %>%
    norm_tss(.) %>% 
    psmelt(.)  %>%
    dplyr::select(c("Sample", "Abundance", taxa_level)) %>%
    pivot_wider(names_from = taxa_level, values_from = "Abundance") %>%
    column_to_rownames(var="Sample") %>%
    mutate(across(where(is.numeric), ~ round(. * scaling_factor))) # Multiply and round

  # Extract row names
  row_names <- rownames(counts_data)
  
  # Convert columns to integer
  counts_data <- as.data.frame(lapply(counts_data, as.integer))
  
  # Reassign row names
  rownames(counts_data) <- row_names
  

  #Reformat data into long format
  counts_data_long <- counts_data %>%
    pivot_longer(
      cols = everything(), 
      names_to = "taxa", 
      values_to = "Count"
    )
  
  #Create a histogram
  p <- ggplot(counts_data_long, aes(x = Count)) +
    geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Histogram of Values", x = "Value", y = "Frequency") +
    stat_bin(
        binwidth = 0.5,
        aes(label = after_stat(count)),
        geom = "text",
        vjust = -0.5, # Adjust the vertical position of the text
        color = "black",
        size = 3
    ) +
    theme_bw()
  
  print(p)
  
  
  # Return metadata and counts data as a named list
  return(list(meta_data = meta_data, counts_data = counts_data))
  
}


### Use FindTopicsNumber() 
#Calculates different metrics to estimate the most preferable number of topics for LDA model.
#CaoJuan2009: https://www.sciencedirect.com/science/article/pii/S092523120800372X
#Arun2010: https://link.springer.com/chapter/10.1007/978-3-642-13657-3_43

#topics <- seq(from = 2, to = 50, by = 2)
#RunFindTopicsNumber(counts_data, topics, "Gibbs")

RunFindTopicsNumber <- function(counts_data, topics, method){ 
  result <- FindTopicsNumber(
    counts_data,
    topics = topics,
    metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
    method = "Gibbs",
    control = list(seed = 1234),
    mc.cores = 2L,
    verbose = TRUE
  ) 
  
  FindTopicsNumber_plot(result)
  p <- FindTopicsNumber_plot(result)

  # Combine each metric into a named row of a data frame
  topic_numbers <- result$topics
  
  df_metrics <- data.frame(
    Metric = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
    stringsAsFactors = FALSE
  )
  
  # Extract each metric as a named row
  df_metrics_values <- rbind(
    setNames(as.data.frame(t(result$Griffiths2004)), paste0("k_", topic_numbers)),
    setNames(as.data.frame(t(result$CaoJuan2009)), paste0("k_",topic_numbers)),
    setNames(as.data.frame(t(result$Arun2010)), paste0("k_",topic_numbers)),
    setNames(as.data.frame(t(result$Deveaud2014)), paste0("k_",topic_numbers))
  )
  
  # Combine metric names and values
  df_result <- cbind(df_metrics, df_metrics_values)
  
  #Find optimial topic number where:
  ##Griffiths2004 and Deveaud2014 are maximized
  ##CaoJuan2009 and Arun2010 are minimized
  
  #Transpose df_result to long format
  df_long <- df_result %>%
    pivot_longer(cols = starts_with("k"), names_to = "k", values_to = "score")
  
  df_long %>%
    group_by(Metric) 
  
  # Normalize each metric (min-max scaling)
  df_scaled <- df_long %>%
    group_by(Metric) %>%
    mutate(score_scaled = (score - min(score)) / (max(score) - min(score))) %>%
    ungroup()
  
  # Flip scores for metrics to minimize
  df_scaled <- df_scaled %>%
    mutate(
      score_final = case_when(
        Metric %in% c("CaoJuan2009", "Arun2010") ~ 1 - score_scaled,
        TRUE ~ score_scaled
      )
    )
  
  # Average score per k across all metrics
  df_combined <- df_scaled %>%
    group_by(k) %>%
    dplyr::summarise(
      combined_score = mean(score_final),
      .groups = "drop"
    ) %>%
    arrange(dplyr::desc(combined_score))
  
  # Best topic number overall
  best_k <- df_combined$k[1]
  cat("Best number of topics (overall balance):", best_k, "\n")
  k_number <- as.numeric(gsub("k", "", best_k))
  
  return(list(full_results = result, plot = p, metrics_df = df_combined, best_k_number = k_number))
}




### LDA: Create the model
#Using LDA to create the topic model 
#https://bookdown.org/josephine_lukito/j381m_tutorial/ldatm.html

#result <- create_topic_model(counts_data, 2, 1, "Gibbs")
#model <- result$lda_model

#CREATE THE MODEL
create_topic_model <- function(counts_data, k_value, alpha, method){
  
    lda_model <- topicmodels::LDA(counts_data, k = k_value, method = method, 
                                  control = list(seed = 1234)) #alpha=alpha add this in if you want to control alpha
    print(paste0("Perplexity is: ", topicmodels::perplexity(lda_model, newdata=as.matrix(counts_data))))
    
    b_df <- data.frame(tidytext::tidy(lda_model, matrix = "beta"))
    g_df <- data.frame(tidytext::tidy(lda_model, matrix = "gamma")) %>%
        arrange(document, topic)
    print(paste0("Alpha value is: ", lda_model@alpha))

    return(list(lda_model = lda_model, b_df = b_df, g_df=g_df))
}



### Plot beta
#plot_beta(result, 15)
#Plot beta, or the numbers that are assigned to each word in a topic. If the beta score is higher, the word matters more in that topic.
plot_beta <- function(lda_result, n_top_topics, b_df, fill_colors = NULL) {
  
  # ----- Top terms across topics for barplots
  top_terms <- b_df %>% 
    group_by(topic) %>% 
    slice_max(beta, n = n_top_topics, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      topic = factor(topic),
      term  = gsub("^g__", "", term),
      term  = str_replace_all(term, "\\.", ""),
      term  = str_replace_all(term, "_", " "),
      term  = str_replace_all(term, "__(1|2|3)$", ""),
      term  = str_replace_all(term, " s  ", " "), 
      term  = gsub("Streptococcus oralis subsp dentisani clade 058", "Streptococcus oralis subsp dentisani", term),
      term  = gsub(" bacterium ", " ", term)
      
    )
  
  # ----- Barplots
  barplots <- ggplot(top_terms, aes(
    x = tidytext::reorder_within(term, beta, topic),
    y = beta,
    fill = topic
  )) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ topic, scales = "free", labeller = labeller(topic = ~ paste0("Topic ", .x)), nrow=1) +
    scale_x_reordered() +
    coord_flip() +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(face = "italic", color="black", size=14),
      axis.text.x = element_text(angle=45, color="black", hjust=1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    labs(x = NULL, y = "Beta")
  
  if (!is.null(fill_colors)) {
    barplots <- barplots + scale_fill_manual(values = fill_colors)
  } else {
    barplots <- barplots + scale_fill_viridis_d()
  }
  
  # Return: barplots ggplot + list of pheatmap objects
  return(barplots)
}

# Heatmap of the relab of the top terms in each topic
plot_beta_heatmaps <- function(lda_result, n_topics, n_top_topics, normalized_count_table, b_df, fill_colors = NULL) {
  
  # ----- One heatmap per topic
  build_heatmap_for_topic <- function(topic_no, n_top = 10) {
    
    # Top terms for this topic (cleaned)
    tt <- b_df %>% 
      filter(topic == topic_no) %>%
      slice_max(beta, n = n_top, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        topic = factor(topic),
        term  = gsub("^g__", "", term),
        term  = str_replace_all(term, "\\.", ""),
        term  = str_replace_all(term, "_", " "),
        term  = str_replace_all(term, "__(1|2|3)$", ""),
        term  = str_replace_all(term, " s  ", " "), 
        term  = gsub("Streptococcus oralis subsp dentisani clade 058", "Streptococcus oralis subsp dentisani", term),
        term  = gsub(" bacterium ", " ", term)
      )
    
    # ---- Matrix for heatmap (with duplicate-handling)
    data_mat <- normalized_count_table %>%
      arrange(Type) %>%
      mutate(Type = as.character(Type)) %>%
      arrange(Type) %>%
      select(-Type) %>%
      tibble::column_to_rownames("Sample") %>%
      as.matrix() %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        term = gsub("^g__", "", term),
        term = stringr::str_replace_all(term, "\\.", ""),
        term = stringr::str_replace_all(term, "_", " "),
        term = stringr::str_replace_all(term, "__(1|2|3)$", "")
      ) %>%
      # If cleaning collapses multiple rows to the same 'term',
      # combine them (sum; use mean if you prefer)
      dplyr::group_by(term) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), sum), .groups = "drop") %>%
      # Keep only the topic's top terms
      dplyr::filter(term %in% tt$term) %>%
      tibble::column_to_rownames("term") %>%
      as.matrix()
    
    # ---- Robust order vector: keep only terms that actually exist as rows
    term_order <- unique(tt$term)
    term_order <- term_order[term_order %in% rownames(data_mat)]
    
    if (length(term_order) == 0) {
      stop("None of the top terms for topic ", topic_no, " were found in normalized_count_table after cleaning.")
    }
    
    # Reorder rows to match top-terms order
    data_mat <- data_mat[term_order, , drop = FALSE]
    
    # Column annotations (samples x Type)
    type_annot <- normalized_count_table %>%
      select(Sample, Type) %>%
      tibble::column_to_rownames("Sample")
    
    ann_colors <- list(Type = c(Plaque = "#3185FC", Abscess = "#FF495C"))
    
    # Build heatmap (silent)
    pheatmap(
      log(data_mat),                  # safer than log()
      annotation_col    = type_annot,
      annotation_colors = ann_colors,
      scale             = "none",
      cluster_rows      = FALSE,        # preserve your order
      cluster_cols      = FALSE,
      color             = colorRampPalette(c("#564592", "white", "#D7CF07"))(100),
      show_rownames     = TRUE,
      show_colnames     = FALSE,
      main              = paste("Topic", topic_no),
      silent            = TRUE
    )
  }
  
  # Build heatmaps for topics 1, 2, 3 using n_top_topics
  hm_list <- lapply(1:n_topics, build_heatmap_for_topic, n_top = n_top_topics)
  
  return(hm_list)
  
}

### Heatmap of gamma scores 
#heatmap_gamma(result)

#If you view the topics_wide data frame, you can see that each document has a gamma score for each topic. Some gamma scores are larger than others. This suggests that a document’s content is predominantly in one topic as opposed to another.

heatmap_gamma <- function(lda_results, type_column, g_df){
    
    topics_wide <- g_df %>%
      pivot_wider(names_from = topic,
                  values_from = gamma)

    topics_wide_type <- meta_data%>% select(type_column) %>% 
      rownames_to_column(var="document") %>% 
      merge.data.frame(topics_wide, by="document") %>%
      dplyr::arrange(type_column) %>%
       mutate(!!type_column := as.character(!!sym(type_column))) %>%  # Convert to character to ensure alphabetical sorting
          dplyr::arrange(!!sym(type_column)) 

    #Set annotation colors
    type_annot <- meta_data %>%
          select(c(type_column))
    ann_colors <- list(Type = c(Plaque = "#3185FC", Abscess = "#FF495C"))
    
    data <- topics_wide_type %>% 
      select(-type_column) %>% 
      remove_rownames() %>% column_to_rownames(var="document") %>% t()
    
    # Create the heatmap
    pheatmap(
          log(data),
          annotation_col = type_annot,  # Add column annotations
          scale = "none",               # Scale rows (optional)
          cluster_rows = FALSE,         # Cluster rows
          cluster_cols = FALSE,         # Cluster columns
          color = colorRampPalette(c( "#564592", "white", "#D7CF07"))(100),  # Custom color palette
          show_rownames = TRUE,        # Show row names
          show_colnames = TRUE,         # Show column names
          annotation_colors = ann_colors,
          gaps_row =  1:(nrow(data) - 1)
        )
}


### UMAP of gamma scores
#type_colors <- c("Abscess" = "red", "Plaque" = "blue")
#plot_gamma_umap(results, "Type", type_colors)


plot_gamma_umap <- function(lda_results, type_column, type_colors, g_df ){
  # Wide gamma matrix
  topics_wide <- g_df %>%
    tidyr::pivot_wider(
      names_from  = topic,
      values_from = gamma
    )
  
  # Join metadata (assumes meta_data has rownames = document IDs)
  res_with_annotations <- meta_data %>%
    tibble::rownames_to_column(var = "document") %>%
    dplyr::select(document, {{ type_column }}) %>%
    merge.data.frame(topics_wide, by = "document")
  
  
  # UMAP config
  custom.config <- umap::umap.defaults
  custom.config$random_state <- 1234
  
  # Run UMAP on topics only
  umap_result <- res_with_annotations %>%
    dplyr::select(-document, -{{ type_column }}) %>%
    umap::umap(config = custom.config)

  # Build a data.frame for ggplot
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$document <- res_with_annotations$document
  umap_df$Type     <- res_with_annotations %>% dplyr::pull({{ type_column }})
  
  # ggplot UMAP
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Type)) +
    geom_point(size = 2, alpha = 0.9) +
    labs(
      title = "UMAP Visualization",
      x     = "UMAP 1",
      y     = "UMAP 2"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(type_colors)) {
    umap_plot <- umap_plot + scale_color_manual(values = type_colors)
  }
  
  # For silhouette, use the same matrix used for UMAP
  pca_matrix <- umap_result$layout
  
  # Use the same Type as in umap_df to stay aligned
  type_factor <- factor(umap_df$Type)
  dist_matrix <- dist(pca_matrix)
  cluster_labels <- as.integer(type_factor)
  
  sil <- cluster::silhouette(cluster_labels, dist_matrix)
  
  sil_df <- as.data.frame(sil[, 1:3])
  sil_df$Type <- type_factor
  colnames(sil_df) <- c("Cluster", "Neighbor", "Silhouette", "Type")
  
  sil_plot <- ggplot(sil_df, aes(x = Type, y = Silhouette, fill = Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    stat_summary(
      fun = median, geom = "text",
      aes(label = round(..y.., 2)),
      vjust = .9, color = "white",
      fontface = "bold", size = 4.5
    ) +
    labs(
      title = "Silhouette Scores by Type",
      x     = "Type",
      y     = "Silhouette Score"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(type_colors)) {
    sil_plot <- sil_plot + scale_fill_manual(values = type_colors)
  }
  
  return(list(umap = umap_plot, sil_plot = sil_plot))
}

###Membership of topic by Type

topic_membership <- function(lda_result, type_column, colors, g_df){
  
    topics_wide <- g_df %>%
      tidyr::pivot_wider(names_from = topic, values_from = gamma)  # columns "1","2",...
    
    topics_long_type <- meta_data %>%
      dplyr::select(all_of(type_column)) %>%
      tibble::rownames_to_column(var = "document") %>%
      dplyr::left_join(topics_wide, by = "document") %>%
      tidyr::pivot_longer(
        cols = -c(document, all_of(type_column)),
        names_to = "topic",
        values_to = "gamma"
      ) %>%
      mutate(topic = factor(topic, levels = sort(unique(topic))))

    # Run Wilcoxon test (default) or t.test
    pvals <- topics_long_type %>%
      dplyr::group_by(topic) %>%
      dplyr::summarise(
        p_value = wilcox.test(gamma ~ Type)$p.value,
        .groups = "drop"
      ) %>%
      mutate(p_adj = p.adjust(p_value, method = "fdr"))  # adjust for multiple testing
    
    
    plot <- ggplot(topics_long_type, aes(x = Type, y = gamma, fill = Type)) +
      geom_boxplot() +
      scale_fill_manual(values = colors) +
      facet_wrap(~ topic, labeller = labeller(topic = ~ paste0("Topic ", .x)), nrow=1) +
      labs(x = "Type", y = expression(gamma)) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) 
    
    return(list(plot = plot))
}




### Heatmap of rel abundance in original data of top taxa
heatmap_by_topic <- function(b_df = result$b_df, normalized_count_table = ps_fs_genus_csv, n_top = 15, topic_no = 1){
  
  # Get top terms ordered by beta (desc)
  top_terms <- b_df %>% 
    dplyr::filter(topic == topic_no) %>%
    dplyr::arrange(desc(beta)) %>%      # <-- order by beta
    dplyr::slice_head(n = n_top) %>%   # safer than top_n()
    dplyr::ungroup() %>%
    mutate(
      topic = factor(topic))

  
  data <- normalized_count_table %>%
    dplyr::arrange(Type) %>%
    mutate(Type = as.character(Type)) %>%    # <- fix `:=` here
    dplyr::arrange(Type) %>% 
    dplyr::select(-Type) %>% 
    column_to_rownames(var = "Sample") %>% 
    as.matrix() %>% 
    t() %>% 
    as.data.frame() %>%
    rownames_to_column(var = "term") %>%
    dplyr::filter(term %in% top_terms$term) %>%
    # bring in beta and sort by it
    dplyr::left_join(top_terms %>% dplyr::select(term, beta), by = "term") %>%
    dplyr::arrange(desc(beta)) %>%
    dplyr::select(-beta) 
  
data <- data %>% mutate(
  term  = gsub("^g__", "", term),
  term  = str_replace_all(term, "\\.", ""),
  term  = str_replace_all(term, "_", " "),
  term  = str_replace_all(term, "__(1|2|3)$", ""),
  term  = str_replace_all(term, " s  ", " "),
  term  = gsub("Streptococcus oralis subsp dentisani clade 058", "Streptococcus oralis", term),
  term  = gsub(" bacterium ", " ", term)
  )%>%
  column_to_rownames(var = "term")
  
  type_annot <- ps_fs_genus_csv %>%
    select(c(Type, Sample)) %>% column_to_rownames(var="Sample")
  
  ann_colors <- list(Type = c(Plaque = "#3185FC", Abscess = "#FF495C"))
  
  italic_labels <- parse(text = paste0("italic('", rownames(data), "')"))
  
  # Create the heatmap
  p <- pheatmap(log(data),
                main = paste0("Topic ", topic_no),
                annotation_col = type_annot,  # Add column annotations
                scale = "none",               # Scale rows (optional)
                fontsize_row = 12,
                cluster_rows = FALSE,         # Cluster rows
                cluster_cols = FALSE,         # Cluster columns
                color = colorRampPalette(c( "#564592", "white", "#D7CF07"))(100),  # Custom color palette
                show_rownames = TRUE,        # Show row names
                show_colnames = FALSE,         # Show column names
                labels_row = italic_labels,       
                annotation_colors = ann_colors 
  )
  return(list(heatmap = p, data =data)) 
}

