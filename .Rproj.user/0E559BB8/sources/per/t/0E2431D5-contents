 
library("pacman")
pacman::p_load("tidyverse")

set.seed(1234)
date = Sys.Date()
 

#Make output directory for data saved from functions
# Specify the directory path
dir_path <- "saved_analysis_files"

# Check if the directory exists
if (!dir.exists(dir_path)) {
  # Create the directory if it doesn't exist
  dir.create(dir_path, recursive = TRUE)  # 'recursive' allows creating parent directories if they don't exist
  message("Directory created: ", dir_path)
} else {
  message("Directory already exists: ", dir_path)
}



# Data Processing

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


## Pediatric: Filter for unmatched samples

#This function is used inside the function filter_phyloseq
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


## Filter for unmatched samples
#This function is used inside the function filter_phyloseq
filter_unmatched_samples <- function(phyloseq_obj) {
  # Extract the sample data
  sample_data_df <- sample_data(phyloseq_obj)
  sample_data_df <- sample_data_df %>% as.matrix() %>% as.data.frame()
  
  # Ensure the necessary columns exist
  required_columns <- c("Sample", "Type")
  if (!all(required_columns %in% colnames(sample_data_df))) {
    stop("The sample data must contain 'Sample' and 'Type' columns.")
  }
  
  # Count the number of unique 'Type' entries for each individual
  type_counts <- sample_data_df %>%
    dplyr::select(c("Sample", "Type")) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarize(num_types = n_distinct(Type), .groups = "drop") 

    # Filter for samples with only one type
  single_type_samples <- type_counts %>%
    filter(num_types == 1) %>%
    pull(Sample)
  
  print(single_type_samples)
  
  filtered_phyloseq_obj <- prune_samples(!(phyloseq_obj@sam_data$Sample) %in% single_type_samples, phyloseq_obj)

  # Return the list of sample IDs
  return(filtered_phyloseq_obj)

}




## Filter phyloseq function
#Example:
#filter_phyloseq(phy_obj_oscc, "OSCC", Contam_g, Contam_f, Contam_s)

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
  
  print(table(phylo_obj_f@sam_data$Sample))
  
  #Remove contaminates
  phylo_obj_f <- subset_taxa(phylo_obj_f, !(Genus %in% Contam_g)) 
  phylo_obj_f <- subset_taxa(phylo_obj_f, !(Family %in% Contam_f)) 
  phylo_obj_f <- subset_taxa(phylo_obj_f, !(Species %in% Contam_s))
    
  #Output
  assign(paste0(phylo_name, "_f"), phylo_obj_f, envir = .GlobalEnv)
}





## Process Phyloseq
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
  
  # Step 7: Dynamically name and assign outputs
  assign(paste0(phylo_name, "_t_n"), phylo_obj_fs_t_n, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_pglom_df"), phylo_obj_fs_t_n_pglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_gglom_df"), phylo_obj_fs_t_n_gglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_sglom_df"), phylo_obj_fs_t_n_sglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_species_csv"), species_csv, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_genus_csv"), genus_csv, envir = .GlobalEnv)
}


#Basic Taxa Plots

## Plot all taxa
#Example: 
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
      filter(!!sym(taxa_level) %in% taxa_filter)
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
      strip.text = element_text(face = "bold")
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
        strip.text = element_text(face = "bold")
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
          strip.text = element_text(face = "bold")
        )
  } else {
      stop("Invalid plot type. Please specify 'box', 'violin', 'density' or 'dot'.")
  }
  print(plot)
}



# Oridinations 
## Plot PCA

#Example (to iterate through a list of ranks):
#resolution <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

#plots <- lapply(resolution, function(rank) {
#  plot_pca(phyloseq_obj = merged_phylo_obj_f, 
#           rank_transformation = rank, 
#           variable = "Study", 
#           colors_list = colors_study)})

#combined_plot <- wrap_plots(plots, ncol = 2) & theme(legend.position = "bottom") 

# Function to create plots
plot_pca <- function(phyloseq_obj, rank_transformation, variable, colors_list=NULL) {
  # Transform and calculate distance
  phylo_trans <- phyloseq_obj %>% tax_fix() %>% tax_transform(rank = rank_transformation, trans = "identity")
  dist_matrix <- dist_calc(phylo_trans, "euclidean")
  ord_res <- ord_calc(dist_matrix, "PCA")
  
  # Plot
  p <- ord_plot(ord_res, axes = c(1, 2), fill = variable, shape = variable, alpha = 0.8, size = 2) +
    ggtitle(label = paste0(
      rank_transformation)) +
    scale_shape_girafe_filled() +
    ggplot2::stat_ellipse(aes(color = !!sym(variable) ) )+
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = .5),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.text = element_text(size = 10, vjust = 0.5, hjust = 1)
    )
  
    # Add color scales if a custom colors_list is provided
  if (!is.null(colors_list)) {
    p <- p + scale_fill_manual(values = colors_list) +
      scale_color_manual(values = colors_list)
  }
  return(p)
}



## Plot PCoA

#Example (to iterate through a list of ranks):
#resolution <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
#plots <- lapply(resolution, function(rank) {
#  plot_PCoA(phyloseq_obj = merged_phylo_obj_f,
#    rank_transformation = rank,
#    trans_type = "identity",        
#    dist_cal_type = "bray",   
#    ord_calc_method = "NMDS",
#    variable = "Study", 
#    colors_list = colors_study)})
    
#combined_plot <- wrap_plots(plots, ncol = 2) & theme(legend.position = "bottom") 

# Function to create plots
plot_PCoA <- function(phyloseq_obj, rank_transformation, trans_type, dist_cal_type, ord_calc_method, variable, colors_list=NULL) {
  # Transform and calculate distance
  phylo_trans <- phyloseq_obj %>% tax_fix() %>%
    tax_transform(rank = rank_transformation, trans = trans_type)
  dist_matrix <- dist_calc(phylo_trans, dist_cal_type)
  ord_res <- ord_calc(dist_matrix, ord_calc_method)
  

  # Plot
  p <- ord_plot(ord_res, axes = c(1, 2), fill = variable, shape = variable, alpha = 0.8, size = 2, plot_taxa = 1:5, size = 2) +
    ggtitle(label = paste0(
      rank_transformation)) +
    scale_shape_girafe_filled() +
    ggplot2::stat_ellipse(aes(color =!!sym(variable) )) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = .5),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.text = element_text(size = 10, vjust = 0.5, hjust = 1)
    )
  
    # Add color scales if a custom colors_list is provided
  if (!is.null(colors_list)) {
    p <- p + scale_fill_manual(values = colors_list) +
      scale_color_manual(values = colors_list)
  }
  return(p)

}



# DA Analyses

## Maaslin2

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

run_Maaslin2 <- function(ps_obj, taxa_level, group, analysis_method, normalization, transform, plot_colors, plot_type, qval_threshold, plot_heights){
  
  # Create a unique hash for the parameters
  param_hash <- digest(list(ps_obj, taxa_level, group, analysis_method, normalization, transform, plot_colors, plot_type, qval_threshold))
  # File path for the analysis results
  result_file <- file.path("saved_analysis_files/", paste0("maaslin2_result_", param_hash, ".rds"))

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
  
  maaslin2_results <- Maaslin2(counts_input , metadata_input, 
                     fixed_effects = group, 
                     plot_heatmap = F, 
                     plot_scatter = F,
                     output = paste0(date, "maaslin2_output"), 
                     analysis_method = analysis_method,
                     normalization = normalization, 
                     transform  = transform)
      # Save the results
      saveRDS(maaslin2_results, result_file)
    }
  
  # Assign the results df to the global environment to output
  assign(paste0("maaslin2_results_", taxa_level),maaslin2_results$results, envir = .GlobalEnv)

  #PLOT EVERYTHING
  
    maaslin_res_filt <- maaslin2_results$results %>% as.data.frame() %>% filter(qval <= qval_threshold) %>%
        dplyr::mutate(Enrichment = ifelse(coef > 0, "Plaque", "Abscess"))
  
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
                axis.text.y = element_text(face="italic"))
    
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
          axis.text.y = element_text(hjust = 1)) 
      

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
      strip.text = element_text(face = "bold")
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
        strip.text = element_text(face = "bold")
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
          strip.text = element_text(face = "bold")
        )
  } else {
      stop("Invalid plot type. Please specify 'box', 'violin', 'density' or 'dot'.")
    }
    
  
    #Plot EVERYTHING and return
    layout <- (bar_plot | dot_plot) / plot + 
    plot_layout(heights = plot_heights) 
    print(layout)
    
  }


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

iterate_maaslin2 <- function(ps_obj, iterative_methods, resolutions, group, qval_threshold, plot_colors, percentage) {
  
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
  result_file <- file.path("saved_analysis_files/", paste0("iterative_maaslin2_result_", param_hash, ".rds"))

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
            output = "maaslin_iterations",
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
    print(p)
    
}







## ANCOM

### Run ANCOM
#Example:

#runANCOM(merged_phylo_obj_fs, 
#         tax_level = "Genus", 
#         group = "condition", 
#         Log2FC_cutoff = 1, 
#         name_of_saved_results = "ancombc_results_genus", 
#         plot_type = "dot", 
#         plot_heights = c(1, 8))

runANCOM <- function(ps_obj, tax_level, group, Log2FC_cutoff=NULL, name_of_saved_results=NULL, plot_type, plot_heights, plot_colors){
  
    # Create a unique hash for the parameters
  param_hash <- digest(list(ps_obj, tax_level, group, Log2FC_cutoff=NULL))
  # File path for the analysis results
  result_file <- file.path("saved_analysis_files/", paste0("ancombc2_result_", param_hash, ".rds"))

  if (file.exists(result_file)) {
      message("Analysis already run. Loading results...")
      ancombc2_results <- readRDS(result_file)
    } else {
      message("Running analysis...")

  # Run ANCOMBC2
  ancombc2_results <- ancombc2(
    data = ps_obj,
    assay_name = "counts",   # Default name for the abundance matrix in phyloseq
    tax_level = tax_level,   # Taxonomic level
    fix_formula = group,     # Fixed effect formula (e.g., "condition")
    p_adj_method = "holm",   # Adjust p-values (e.g., Holm method)
    group = group,           # Grouping variable
    struc_zero = TRUE,       # Detect structural zeros
    neg_lb = TRUE,           # Detect sampling zeros
    alpha = 0.05,            # Significance level
    n_cl = 1                 # Number of cores for parallelization
  )
  
  # Save the results
  saveRDS(ancombc2_results, result_file)
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
        print(paste("Column found:", dff_col_name))
      } else {
        print("No matching column found.")
      }
      
      #LFC column
      lfc_col_name <- grep(paste0("^lfc_", group), colnames(df_fig), value = TRUE)
      # If you want to check and store it
      if (length(lfc_col_name) > 0) {
        # Save the first matching column name (if multiple matches)
        lfc_col_name <- lfc_col_name[1]
        print(paste("Column found:", lfc_col_name))
      } else {
        print("No matching column found.")
      }
      
      #SS condition column
      ss_col_name <- grep(paste0("^passed_ss_", group), colnames(df_fig), value = TRUE)
      # If you want to check and store it
      if (length(ss_col_name) > 0) {
        # Save the first matching column name (if multiple matches)
        ss_col_name <- ss_col_name[1]
        print(paste("Column found:", ss_col_name))
      } else {
        print("No matching column found.")
      }
      
      #SE condition column
      se_col_name <- grep(paste0("^se_", group), colnames(df_fig), value = TRUE)
      # If you want to check and store it
      if (length(se_col_name) > 0) {
        # Save the first matching column name (if multiple matches)
        se_col_name <- se_col_name[1]
        print(paste("Column found:", se_col_name))
      } else {
        print("No matching column found.")
      }
      
      #P-val condition column
      p_col_name <- grep(paste0("^p_", group), colnames(df_fig), value = TRUE)
      # If you want to check and store it
      if (length(p_col_name) > 0) {
        # Save the first matching column name (if multiple matches)
        p_col_name <- p_col_name[1]
        print(paste("Column found:", p_col_name))
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
    print(ancombc2_results$res)
    stop("There are no results to display.")
  }
  
  print(df_fig)

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
      axis.text.y = element_text(hjust = 1, color = df_fig$color)) 

  
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
        axis.text.y = element_text(hjust = 1, color = df_fig$color)) +
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
      strip.text = element_text(face = "bold")
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
        strip.text = element_text(face = "bold")
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
          strip.text = element_text(face = "bold")
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
            strip.text = element_text(face = "bold")
          )
  

  } else {
    stop("Invalid plot type. Please specify 'box', 'violin', 'density' or 'dot'.")
  }
  

  #Plot EVERYTHING and return
  layout <- (bar_plot | dot_plot) / plot + 
  plot_layout(heights = plot_heights) 
  layout
}


## ALDEX2

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
  param_hash <- digest(list(ps_obj, resolutions, group))
  # File path for the analysis results
  result_file <- file.path("saved_analysis_files/", paste0("iterative_aldex2_result_", param_hash, ".rds"))

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
  
  print(summary_table)
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
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
   # scale_fill_manual(values = plot_colors)
  }
    print(p)
  
}



### Run Aldex


run_aldex2 <- function(ps_obj, group, tax_rank, method, transform, normalization, plot_colors){

  aldex2_results <- microbiomeMarker::run_aldex(
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
  
  aldex2_res <- aldex2_results@marker_table %>% as.matrix() %>% as.data.frame()
  assign(paste0("aldex2_res_", tax_rank), aldex2_res, envir = .GlobalEnv)

  #The output of Aldex2 here is kind of confusing, so I am going to make a new column to specify the enrichment
  df_fig <- aldex2_res %>% 
      mutate(ef_aldex = as.numeric(ef_aldex)) %>%
      mutate(padj = as.numeric(padj)) %>%
      dplyr::mutate(Enrichment = ifelse(ef_aldex > 0, "Plaque", "Abscess"))

  df_fig$feature = factor(df_fig$feature, levels = df_fig$feature)
  df_fig$Enrichment = factor(df_fig$Enrichment, 
                           levels = c("Plaque", "Abscess"))

    #DOT PLOT
  dot_plot <- df_fig %>%
    ggplot(aes(x = abs(ef_aldex), y = feature)) +
    geom_point(aes(size = -log10(padj), color = Enrichment)) + # Map size and color
    scale_color_manual(values = plot_colors) +
    scale_size(range = c(1, 5), name = "-log10(P-value)") + # Adjust dot size and legend label
    labs(title = NULL, x="EF Aldex", y = NULL, color = "Enrichment") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
      plot.title = element_text(hjust = 0.5))  # Center the title

  
  #Bar plot
  bar_plot <- df_fig %>%
    ggplot(aes(x = feature, y = abs(ef_aldex), fill = Enrichment)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4))+
    labs(x = NULL, y = "EF Aldex", 
         title = NULL) + 
    scale_fill_manual(values = plot_colors) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank()) +
    coord_flip()
  
  layout <- (bar_plot | dot_plot) 
  layout
  
}


## Combine DA Methods

combine_DA <- function(maaslin2_results, ancombc2_results, aldex2_results, group, tax_level, plot_heights, plot_colors){

    # Prepare the two dataframes for merging
  #MAASLIN2
    maaslin2_results_clean <- maaslin2_results %>%
      dplyr::rename(taxon = feature) %>% # Rename 'feature' to 'taxon' for consistency
      select(taxon, Maaslin2_coef = coef, Maaslin2_pval = pval, Maaslin2_value = value)
    maaslin2_results_clean$taxon <- gsub("s__","SS", maaslin2_results_clean$taxon)
    maaslin2_results_clean$taxon <- gsub("[[:punct:]]","", maaslin2_results_clean$taxon)

    #ANCOMBC2
    dff_col_name <- grep(paste0("^diff_", group), colnames(ancombc2_results), value = TRUE)
    lfc_col_name <- grep(paste0("^lfc_", group), colnames(ancombc2_results), value = TRUE)
    ss_col_name <- grep(paste0("^passed_ss_", group), colnames(ancombc2_results), value = TRUE)
    p_col_name <- grep(paste0("^p_", group), colnames(ancombc2_results), value = TRUE)
          
    ancombc2_results_clean <- ancombc2_results %>%
            select(taxon, ANCOMBC2_DFF = dff_col_name, ANCOMBC2_LFC = lfc_col_name, ANCOMBC2_p_val =p_col_name, ANCOMBC2_passed_SS = ss_col_name) %>%  # Keep relevant columns
            filter(ANCOMBC2_DFF == TRUE) 
    ancombc2_results_clean$taxon <- gsub("s__","SS", ancombc2_results_clean$taxon)
    ancombc2_results_clean$taxon <- gsub("[[:punct:]]","", ancombc2_results_clean$taxon)


    #ALDEX2
      #The output of Aldex2 here is kind of confusing, so I am going to make a new column to specify the enrichment
    aldex2_res_clean <- aldex2_results%>% 
          mutate(ef_aldex = as.numeric(ef_aldex)) %>%
          mutate(padj = as.numeric(padj)) %>%
          mutate(feature = as.character(feature)) %>%
          dplyr::mutate(Enrichment = ifelse(ef_aldex > 0, "Plaque", "Abscess")) %>%
        dplyr::rename(taxon = feature) %>% # Rename 'feature' to 'taxon' for consistency
        select(taxon, Aldex2_Enrichment = Enrichment, Aldex2_EF = ef_aldex, Aldex2_padj = padj)
    aldex2_res_clean$taxon <- gsub("s__","SS", aldex2_res_clean$taxon)
    aldex2_res_clean$taxon <- gsub("[[:punct:]]", "", aldex2_res_clean$taxon)


    # Perform a full join on the 'taxon' column
    DA_results_df <- dplyr::full_join(maaslin2_results_clean, ancombc2_results_clean, aldex2_res_clean, by = "taxon") %>%
      dplyr::full_join(aldex2_res_clean, by = "taxon") %>%
      mutate(confidence = case_when(
        rowSums(is.na(select(., Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) == 2 ~ "Low",
        rowSums(is.na(select(., Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) == 1 ~ "Medium",
        rowSums(is.na(select(., Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) == 0 ~ "High"
      )) %>% 
      rowwise() %>%
      mutate(all_positive_or_negative = all(na.omit(c_across(c(Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) > 0) | 
                                        all(na.omit(c_across(c(Maaslin2_coef, ANCOMBC2_LFC, Aldex2_EF))) < 0)) %>%
      ungroup()
    
        if (any(DA_results_df$all_positive_or_negative == FALSE, na.rm = TRUE)) {
        warning("Some values in 'all_positive_or_negative' are FALSE!")}
      
    #Remove the "g" and "s" that are followed by an uppercase (for genus and species)
    DA_results_df$taxon <- gsub("g(?=[A-Z])", "", DA_results_df$taxon, perl = TRUE)
    DA_results_df$taxon <- gsub("SS", " ", DA_results_df$taxon, perl = TRUE)
    
    # View the resulting dataframe
    print(DA_results_df)
    return(DA_results_df)
    
    assign(paste0("DA_results_", tax_level), DA_results_df, envir = .GlobalEnv)
    
    #Create a barplot of Maaslin2 results based on confidence
    #MEDIUM
    df_fig_med <- DA_results_df %>%
      filter(confidence == "Medium") %>%
      mutate(Enrichment =ifelse(Maaslin2_coef > 0, "Plaque", "Abscess")) %>%
      arrange(Enrichment, abs(Maaslin2_coef)) %>% 
      mutate(feature_order = factor(taxon, levels = unique(taxon))) 
    
    bar_plot_med <- df_fig_med %>%
      ggplot(aes(x = feature_order, y = abs(Maaslin2_coef), fill = Enrichment)) +
        geom_bar(stat = "identity", color="black") +
        coord_flip() +
        labs(title = "Medium Confidence", x="", y = "Maaslin2 coef") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 1),
          axis.text.y = element_text(face="italic"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold")
        ) +
        scale_y_discrete(expand = c(0, 0.5)) + # change additive expansion from default 0.6 to 0.5
        scale_fill_manual(values = plot_colors) +
        theme(strip.text = element_text(size = 12), strip.background = element_rect(fill = "lightgrey"))
  
    #HIGH
    df_fig_high <- DA_results_df %>%
      filter(confidence == "High") %>%
      mutate(Enrichment =ifelse(Maaslin2_coef > 0, "Plaque", "Abscess")) %>%
      arrange(Enrichment, abs(Maaslin2_coef)) %>% 
      mutate(feature_order = factor(taxon, levels = unique(taxon))) 
    
    bar_plot_high <- df_fig_high %>%
      ggplot(aes(x = feature_order, y = abs(Maaslin2_coef), fill = Enrichment)) +
        geom_bar(stat = "identity", color="black") +
        coord_flip() +
        labs(title = "High Confidence", x="", y = "Maaslin2 coef") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 1),
          axis.text.y = element_text(face="italic"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(face = "bold")
        ) +
        scale_y_discrete(expand = c(0, 0.5)) + # change additive expansion from default 0.6 to 0.5
        scale_fill_manual(values = plot_colors) +
        theme(strip.text = element_text(size = 12), strip.background = element_rect(fill = "lightgrey"))
    
    layout <- bar_plot_high / bar_plot_med + 
    plot_layout(heights = plot_heights) 
    print(layout)
}



#LDA model function

### Iteratively find thresholds
#threshold_iterations <-  seq(.5, 1, by = .05)
#result_table <- find_thresholds_iteratively(ps_obj_ped_fs, "Genus", threshold_iterations)

find_thresholds_iteratively <- function(psobj, rank, required_percent_taxa_seq) {
  
  # Initialize an empty data frame to store results
  results <- data.frame(
    required_percent_taxa = numeric(),
    number_of_numeric_columns = integer(),
    threshold = numeric(),
    stringsAsFactors = FALSE
  )

  # Agglomerate at rank
  core_ps <- tax_glom(psobj, taxrank = rank)
  core_ps <- norm_tss(core_ps)
    
  # Extract otu matrix and other data (some of which aren't necessary for below, but I'm keeping for future use)
  otu_mat <- as.matrix(otu_table(core_ps))
  total_ASVs <- length(rowSums(otu_mat))
  non_zero_ASVs <- sum(rowSums(otu_mat) > 0)
    
  # Find the minimum non-zero value in the entire OTU table
 min_value <- min(as.vector(otu_mat)[as.vector(otu_mat) > 0], na.rm = TRUE)
    
    # Prepare the taxa matrix for filtering
  tax_mat <- psmelt(core_ps) %>%
      dplyr::select(c("Sample", "Abundance", rank)) %>%
      pivot_wider(names_from = rank, values_from = "Abundance")
  
  # Loop through each value of required_percent_taxa
  for (required_percent_taxa in required_percent_taxa_seq) {

    # Find columns that meet the criteria
    filtered_df <- tax_mat[, colMeans(tax_mat >= min_value, na.rm = TRUE) >= required_percent_taxa]
    
    # Ensure the dataframe has stuff in it
    if (ncol(filtered_df) > 1) {
      
        numeric_cols <- filtered_df[, sapply(filtered_df, is.numeric), drop = FALSE]
        # Count the number of numeric columns
        num_numeric_cols <- ncol(numeric_cols)
        
        # Find the new minimum non-zero value across all numeric columns
        threshold <- min(numeric_cols[numeric_cols > 0], na.rm = TRUE)
    } else {
        # If no numeric columns, set values accordingly
        num_numeric_cols <- 0
        threshold <- NA
    }
    
    # Add the result to the results data frame
    results <- rbind(results, data.frame(
      required_percent_taxa = required_percent_taxa,
      number_of_numeric_columns = num_numeric_cols,
      threshold = threshold,
      stringsAsFactors = FALSE
    ))
  }
  
  print(results)
  # Return the final results table
  return(results)
}



### Find filtering threshold - set value

#threshold <- find_threshold(ps_obj_ped_fs, "Genus", .95)


find_threshold <- function(psobj, rank, required_percent_taxa){

        #Agglomerate at rank
        core_ps <- tax_glom(psobj, taxrank = rank)
        core_ps <- norm_tss(core_ps)

        #Extract otu matrix and other data (some of which aren't necessary for below, but I'm keeping for future use)
        otu_mat <- as.matrix(otu_table(core_ps))
        total_ASVs <- length(rowSums(otu_mat))
        non_zero_ASVs <- sum(rowSums(otu_mat) > 0)
        
        #Find the minimum non-zero value in the entire OTU table, which is what we will pass to the filtering script
        min_value <- min(as.vector(otu_mat)[as.vector(otu_mat) > 0], na.rm = TRUE)
        
        #Prepare the taxa matrix for filtering
        tax_mat <- psmelt(core_ps) %>%
            dplyr::select(c("Sample", "Abundance", rank)) %>%
            pivot_wider(names_from = rank, values_from = "Abundance")
        
        # Find columns that meet the criteria: At least x% (defined in the function by user) of the samples have these taxa above a minimum threshold (calculated above)
        filtered_df <- tax_mat[, colMeans(tax_mat >= min_value, na.rm = TRUE) >= required_percent_taxa]
        
        print(filtered_df)
      if (ncol(filtered_df) > 1) {
          # Select only numeric columns
          numeric_cols <- filtered_df[sapply(filtered_df, is.numeric)]
          
          # Find the new minimum non-zero value across all numeric columns
          threshold <- min(numeric_cols[numeric_cols > 0], na.rm = TRUE)
      } else{
        threshold = NA
        print("No taxa retained at this filtering threshold!")
      }
        # Return that threshold for later use
        return(threshold)
}

 



### Prep and binarize input
#This works for all topic modeling algorithms as a pre-processing step
#results_threshold <- prep_data_binarize(ps_obj_ped_fs, "Genus", threshold)
#meta_data <- results_threshold$meta_data
#counts_data <- results_threshold$counts_data

#Prepare data from phyloseq object and binarize at a certain threshold of relative abundance
prep_data_binarize <- function(phy_obj, taxa_level, binarization_threshold, type_column){
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
    mutate(across(where(is.numeric), ~ ifelse(. >= binarization_threshold, 1, 0))) #Apply binarization threshold

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
#alculates different metrics to estimate the most preferable number of topics for LDA model.
#CaoJuan2009: https://www.sciencedirect.com/science/article/pii/S092523120800372X
#Arun2010: https://link.springer.com/chapter/10.1007/978-3-642-13657-3_43

#topics <- seq(from = 2, to = 50, by = 2)
#RunFindTopicsNumber(counts_data, topics, "Gibbs")

RunFindTopicsNumber <- function(counts_data, topics, method){ 
  result <- FindTopicsNumber(
      counts_data,
      topics = topics,
      metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
      method = method,
      control = list(seed = 1234),
      mc.cores = 2L,
      verbose = TRUE
    ) 
  FindTopicsNumber_plot(result)
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

plot_beta <- function(lda_result, n_top_topics){
    top_terms <- result$b_df %>% 
      group_by(topic) %>% 
      top_n(n_top_topics, beta) %>%
      ungroup() %>%
      mutate(
        term_clean = term %>%
          gsub("^g__|^s__", "", .) %>%      # remove g__ or s__ at the start
          gsub("__[0-9]+$", "", .) %>%       # remove trailing _1, _2, etc.
          gsub("_", " ", .)                 # replace _ with space
      )
    
    top_terms %>% 
      ggplot(aes(
        x = tidytext::reorder_within(term_clean, beta, topic),  # descending order
        y = beta,
        fill = factor(topic)
      )) +
      geom_bar(stat = 'identity', show.legend = FALSE) +
      facet_wrap(~ topic, scales = "free") +
      coord_flip() +
      scale_x_reordered(
        labels = function(x) parse(text = paste0("italic('", x, "')"))
      ) +
      theme_bw(base_size = 12) +
      scale_fill_viridis_d() +
      labs(x = NULL, y = "Beta")
}



### Heatmap of gamma scores 
#heatmap_gamma(result)

#If you view the topics_wide data frame, you can see that each document has a gamma score for each topic. Some gamma scores are larger than others. This suggests that a document’s content is predominantly in one topic as opposed to another.

heatmap_gamma <- function(lda_results, type_column){
    
    topics_wide <- lda_results$g_df %>%
      pivot_wider(names_from = topic,
                  values_from = gamma)

    topics_wide_type <- meta_data%>% select(type_column) %>% 
      rownames_to_column(var="document") %>% 
      merge.data.frame(topics_wide, by="document") %>%
      dplyr::arrange(type_column) %>%
       mutate(!!type_column := as.character(!!sym(type_column))) %>%  # Convert to character to ensure alphabetical sorting
          dplyr::arrange(!!sym(type_column)) 

    annotations <- topics_wide_type %>%
          select(c("document", type_column)) %>%
          remove_rownames() %>%
          column_to_rownames(var = "document")
    
    data <- topics_wide_type %>% 
      select(-type_column) %>% 
      remove_rownames() %>% column_to_rownames(var="document") %>% t()
    
    # Create the heatmap
    pheatmap(
          data,
          annotation_col = annotations,  # Add column annotations
          scale = "row",               # Scale rows (optional)
          cluster_rows = TRUE,         # Cluster rows
          cluster_cols = FALSE,         # Cluster columns
          color = colorRampPalette(c("blue", "white", "red"))(50),  # Custom color palette
          show_rownames = TRUE,        # Show row names
          show_colnames = TRUE         # Show column names
        )
}





### UMAP of gamma scores
#type_colors <- c("Abscess" = "red", "Plaque" = "blue")
#plot_gamma_umap(results, "Type", type_colors)


plot_gamma_umap <- function(lda_results, type_column, type_colors ){
    topics_wide <- result$g_df %>%
          pivot_wider(names_from = topic,
                      values_from = gamma) 
    
    res_with_annotations <- meta_data %>% select(type_column) %>% 
          rownames_to_column(var="document") %>% 
          merge.data.frame(topics_wide, by="document") 
    
    #Set UMAP seed
    custom.config <- umap.defaults
    custom.config$random_state <- 1234
    
    # Run UMAP
    umap_result <- res_with_annotations %>% select(-type_column) %>% column_to_rownames(var="document") %>% umap(., config=custom.config)
    
    # Convert the annotations to numeric for coloring
    annotations <- res_with_annotations %>% column_to_rownames(var="document") %>%  pull(type_column)  # Extract "Type" as a vector
    
    # Assign colors to each type
    colors <- type_colors[annotations]  # Match colors to each sample
    
    # Plot UMAP
    plot(umap_result$layout, 
         col = colors,
         pch = 19,           # Point type
         xlab = "UMAP 1", 
         ylab = "UMAP 2", 
         main = "UMAP Visualization")
}


###Membership of topic by Type

topic_membership <- function(lda_result, type_column, colors){
  
      topics_wide <- lda_result$g_df %>%
            pivot_wider(names_from = topic,
                        values_from = gamma)
      
      topics_long_type <- meta_data%>% select(type_column) %>% 
            rownames_to_column(var="document") %>% 
            merge.data.frame(topics_wide, by="document") %>%
        pivot_longer(cols = -c(document, type_column), 
                     names_to = "topic", 
                     values_to = "gamma")

      topics_long_type %>%
        ggplot(aes(type_column, gamma, fill=!!sym(type_column))) +
        geom_boxplot() +
        scale_fill_manual(values= colors) +
        facet_wrap(~ topic) +
        labs(x = "topic", y = expression(gamma)) +
        theme_bw()
}




### Heatmap of rel abundance in original data of top taxa
#This is throwing an error with some ps objects, so I need to troubleshoot

relab_heatmap <- function(lda_results, psobj, rank, type_column, topic_no, n_top_words){

    #Gather a list of all the most important taxa(words) in the topic
    taxa_list <- lda_results$b_df %>%
      filter(topic == topic_no) %>%
      top_n(n_top_words, beta) %>% #takes the words with the top 10 beta scores
      distinct(term) %>%
      pull(term)
    
    #Agglomerate your phyloseq object at the desired level
    core_ps <- tax_glom(psobj, taxrank = rank)
    core_ps <- norm_tss(core_ps)
    
    tax_mat <- psmelt(core_ps) %>%
                dplyr::select(c("Sample", "Abundance", rank)) %>%
                pivot_wider(names_from = rank, values_from = "Abundance") %>%
                select(any_of(taxa_list), Sample)

    tax_mat_type <- meta_data%>% select(type_column) %>% 
          rownames_to_column(var="Sample") %>% 
          merge.data.frame(tax_mat, by="Sample")%>%
          mutate(!!type_column := as.character(!!sym(type_column))) %>%  # Convert to character to ensure alphabetical sorting
          dplyr::arrange(!!sym(type_column)) 
    
        annotations <- tax_mat_type %>%
          select(c("Sample", type_column)) %>%
          remove_rownames() %>%
          column_to_rownames(var = "Sample")
            
    data <- tax_mat_type %>% 
      select(-type_column) %>% 
      remove_rownames() %>% column_to_rownames(var="Sample") %>% t()
    
    if (ncol(data) == 0 || any(!is.finite(range(data, na.rm = TRUE)))) {
      stop("Heatmap input matrix has no finite values. Check taxonomic rank or taxa list.")
    }
    
    pheatmap(log10((100*data) + .00001),
              annotation_col = annotations,  # Add column annotations
              #scale = "row",               # Scale rows (optional)
              cluster_rows = FALSE,         # Cluster rows
              cluster_cols = FALSE,         # Cluster columns
              color = colorRampPalette(c("blue", "white", "red"))(50),  # Custom color palette
              show_rownames = TRUE,        # Show row names
              show_colnames = FALSE,    
             main = paste0("Topic ", topic_no)
    )  
}


