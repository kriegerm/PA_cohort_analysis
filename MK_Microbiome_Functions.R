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


## Filter for unmatched samples - non-Pediatric study
#This function is used inside the function filter_phyloseq
filter_unmatched_samples <- function(phyloseq_obj) {
  
  #Make sure only Plaque and Abscess samples are retained
  phyloseq_obj <- subset_samples(phyloseq_obj, Type %in% c("Abscess", "Plaque"))
  
  # Extract sample data with sample_names as a column
  sample_data_df <- as(sample_data(phyloseq_obj), "data.frame")
  sample_data_df$SampleID <- rownames(sample_data_df)  # These are sample_names()
  
  # Check for required columns
  if (!all(c("Sample", "Type") %in% colnames(sample_data_df))) {
    stop("Sample data must contain 'Sample' and 'Type' columns.")
  }
  
  #Find Samples (i.e., individuals) that only have one Type
  type_counts <- sample_data_df %>%
    dplyr::filter(Type %in% c("Plaque", "Abscess")) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(num_types = n_distinct(Type), .groups = "drop")
  
  # Get the Sample IDs (individuals) to drop
  single_type_sample_values <- type_counts %>%
    dplyr::filter(num_types == 1) %>%
    dplyr::pull(Sample)
  
  #Find the sample_names (SampleID) associated with those individuals
  sample_names_to_remove <- sample_data_df %>%
    dplyr::filter(Sample %in% single_type_sample_values) %>%
    dplyr::pull(SampleID)
  
  # Print for confirmation
  print("Samples being removed (sample_names):")
  print(sample_names_to_remove)
  
  # Prune phyloseq object
  filtered_phyloseq_obj <- prune_samples(!(sample_names(phyloseq_obj) %in% sample_names_to_remove), phyloseq_obj)
  return(filtered_phyloseq_obj)
}



## ---- Filter Phyloseq ----
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
  
  # Step 7: Dynamically name and assign outputs
  assign(paste0(phylo_name, "_t_n"), phylo_obj_fs_t_n, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_pglom_df"), phylo_obj_fs_t_n_pglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_gglom_df"), phylo_obj_fs_t_n_gglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_t_n_sglom_df"), phylo_obj_fs_t_n_sglom_df, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_species_csv"), species_csv, envir = .GlobalEnv)
  assign(paste0(phylo_name, "_genus_csv"), genus_csv, envir = .GlobalEnv)
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
  return(list(plot = plot))
}

# ---- Ordinations ----

## ---- Ordinations with Microviz ----

#Example (to iterate through a list of ranks):
#resolution <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
#plots <- lapply(resolution, function(rank) {
#  plot_pca(phyloseq_obj = merged_phylo_obj_f, 
#           rank_transformation = rank, 
#           variable = "Study", 
#           colors_list = colors_study)})
#combined_plot <- wrap_plots(plots, ncol = 2) & theme(legend.position = "bottom") 

plot_pca_microviz <- function(phyloseq_obj, rank_transformation, variable, colors_list=NULL) {
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

plot_PCoA_microviz <- function(phyloseq_obj, rank_transformation, trans_type, dist_cal_type, ord_calc_method, variable, colors_list=NULL) {
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

## ---- Compare methods ----

#BIG function to compare different ordination methods 
compare_ordination_methods <- function(phyloseq_obj,
                                       rank_transformations = c("Species", "Genus", "Phylum"),
                                       trans_types = c("relab", "clr", "identity", "log"),
                                       dist_cal_types = c("euclidean", "bray", "jaccard"),
                                       ord_calc_methods = c("PCA", "PCoA", "NMDS"),
                                       No_axes_to_correlated_with_LibSize = 2,
                                       variance_threshold = 0.80) { #find the number of axes needed for PCA and PCoA that explain this percent of variance  
  library(phyloseq)
  library(vegan)
  library(cluster)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(digest)
  library(ggplot2)
  
  #Create a parameter hash so that you can just load results if they've already been run in an identical way
  param_hash <- digest(list(
    study = cfg$study,
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
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(rlang)
  
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
    geom_point(aes(fill = !!sym(variable), color = !!sym(variable)), alpha = 0.8, size = 2, shape = 21) +
    ggtitle(paste0(rank_transformation, " - ", ord_calc_method, " (", trans_type, ", ", dist_cal_type, ")")) +
    stat_ellipse(aes(color = !!sym(variable))) +
    theme_bw() +
    labs(x = x_lab, y = y_lab)+
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = .5),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.text = element_text(size = 10, vjust = 0.5, hjust = 1)
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

## ---- Analyze Sil Score ----

# Analyze Silhouette score using "Type" as cluster label and return plot
analyze_type_clustering_on_pca <- function(scores_df,
                                           component_num = 3,
                                           colors_list = NULL) {
  library(ggplot2)
  library(dplyr)
  library(cluster)
  
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
    geom_boxplot(outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                 vjust = -0.5, color = "white", fontface = "bold", size = 3.5) +
    labs(title = "Silhouette Scores by Type",
         x = "Type",
         y = "Silhouette Score") +
    theme_bw() 
  
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
  param_hash <- digest(list(ps_obj, resolutions, group))
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
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
      plot.title = element_text(hjust = 0.5)
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
      panel.grid.minor.y = element_blank()
    )
  
  #Return
  return(list(results = df_fig, bar_plot = bar_plot, dot_plot = dot_plot))
}


## ---- Combine DA Methods----

combine_DA <- function(maaslin2_results, ancombc2_results, aldex2_results, group, tax_level, plot_colors){

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
    
    return(list(results = DA_results_df, bar_plot_high_confidence = bar_plot_high, bar_plot_med_confidence = bar_plot_med))
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
  
  return(list(full_results = result, metrics_df = df_combined, best_k_number = k_number))
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

plot_beta <- function(lda_result, n_top_topics, b_df){
    top_terms <- b_df %>% 
      group_by(topic) %>% 
      top_n(n_top_topics, beta) %>%
      ungroup() 

    top_terms %>% 
      ggplot(aes(
        x = tidytext::reorder_within(term, beta, topic),  # descending order
        y = beta,
        fill = factor(topic)
      )) +
      geom_bar(stat = 'identity', show.legend = FALSE) +
      facet_wrap(~ topic, scales = "free") +
      coord_flip() +
      # scale_x_reordered(
      #   labels = function(x) parse(text = paste0("italic('", x, "')"))
      # ) +
      theme_bw(base_size = 12) +
      scale_fill_viridis_d() +
      labs(x = NULL, y = "Beta")
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


plot_gamma_umap <- function(lda_results, type_column, type_colors, g_df ){
    topics_wide <- g_df %>%
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

topic_membership <- function(lda_result, type_column, colors, g_df){
  
      topics_wide <- g_df %>%
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

relab_heatmap <- function(lda_results, psobj, rank, type_column, topic_no, n_top_words, b_df){

    #Gather a list of all the most important taxa(words) in the topic
    taxa_list <- b_df %>%
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


