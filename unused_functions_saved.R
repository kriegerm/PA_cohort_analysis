
compare_ordination_methods <- function(phyloseq_obj,
                                       rank_transformations = c( "Species", "Genus", "Phylum"),
                                       trans_types = c("relab", "clr", "identity", "log"),
                                       dist_cal_types = c("euclidean", "bray", "jaccard"),
                                       ord_calc_methods = c("PCA", "PCoA", "NMDS"),
                                       component_num = 2,
                                       k_range = 2:10) {
  library(phyloseq)
  library(vegan)
  library(cluster)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(digest)
  
  
  # Create a unique hash for the parameters
  param_hash <- digest(list(
    phyloseq_obj = phyloseq_obj,
    rank_transformations = rank_transformations,
    trans_types = trans_types,
    dist_cal_types = dist_cal_types,
    ord_calc_methods = ord_calc_methods,
    component_num = component_num,
    k_range = k_range
  ))
  
  # Define result path
  dir.create("saved_analysis_files", showWarnings = FALSE)  # Create folder if needed
  result_file <- file.path("saved_analysis_files", paste0("ordination_parameter_sweep_result_", param_hash, ".rds"))
  
  # Check cache
  if (file.exists(result_file)) {
    message("Analysis already run. Loading results...")
    return(readRDS(result_file))
  } else {
    message("Running analysis...")
  }
  
  clr_transform <- function(mat, pseudocount = 1) {
    mat <- mat + pseudocount
    log_mat <- log(mat)
    gm <- rowMeans(log_mat)
    sweep(log_mat, 1, gm)
  }
  
  validate_transformation_distance <- function(trans_type, dist_type, ord_method) {
    if (trans_type == "clr" && dist_type != "euclidean") {
      stop("CLR is only compatible with euclidean distance.")
    }
    return(dist_type)
  }
  
  results_list <- list()
  counter <- 1
  
  for (rank_trans in rank_transformations) {
    for (trans_type in trans_types) {
      for (ord_method in ord_calc_methods) {
        ord_method_upper <- toupper(ord_method)
        
        if (ord_method_upper == "PCA") {
          dist_types_to_use <- "euclidean"
        } else {
          dist_types_to_use <- dist_cal_types
        }
        
        for (dist_type in dist_types_to_use) {
          message(sprintf("Trying: %s | %s | %s | %s", rank_trans, trans_type, dist_type, ord_method))
          try({
            phylo_trans <- tax_glom(phyloseq_obj, taxrank = rank_trans)
            otu_mat <- otu_table(phylo_trans)
            if (taxa_are_rows(phylo_trans)) {
              otu_mat <- t(otu_mat)
            }
            otu_mat <- as.matrix(otu_mat)
            
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
            
            otu_trans <- otu_trans[, apply(otu_trans, 2, var) != 0]
            dist_type_valid <- validate_transformation_distance(trans_type, dist_type, ord_method)
            
            scores_df <- NULL
            var_expl <- NULL
            stress <- NA
            ord_method_upper <- toupper(ord_method)
            if (ord_method_upper == "PCA") {
              ord <- prcomp(otu_trans, center = TRUE, scale. = TRUE)
              scores_df <- as.data.frame(ord$x[, 1:component_num])
              colnames(scores_df) <- paste0("Axis", 1:component_num)
              var_expl <- ord$sdev^2 / sum(ord$sdev^2)
              
            } else if (ord_method_upper == "PCOA") {
              dist_mat <- vegdist(otu_trans, method = dist_type_valid)
              ord <- cmdscale(dist_mat, k = component_num, eig = TRUE)
              scores_df <- as.data.frame(ord$points)
              colnames(scores_df) <- paste0("Axis", 1:component_num)
              var_expl <- ord$eig / sum(ord$eig)
              
            } else if (ord_method_upper == "NMDS") {
              dist_mat <- vegdist(otu_trans, method = dist_type_valid)
              ord <- metaMDS(dist_mat, k = component_num, trymax = 100)
              scores_df <- as.data.frame(ord$points)
              colnames(scores_df) <- paste0("Axis", 1:ncol(scores_df))
              stress <- ord$stress
            }
            
            scores_df$SampleID <- rownames(scores_df)
            
            meta_df <- sample_data(phylo_trans) %>% as.data.frame() 
            scores_df <- dplyr::left_join(scores_df, meta_df, by = "SampleID")
            scores_df$LibrarySize <- sample_sums(phyloseq_obj)[scores_df$SampleID]
            
            
            
            rho <- as.numeric(suppressWarnings(cor.test(scores_df$LibrarySize, scores_df$Axis1, method = "spearman")$estimate))
            
            
            var1 <- var2 <- NA
            if (!is.null(var_expl)) {
              var1 <- var_expl[1] * 100
              var2 <- var_expl[2] * 100
            }
            
            pca_mat <- scores_df[, paste0("Axis", 1:component_num)]
            sil_scores <- map_dbl(k_range, function(k) {
              km <- kmeans(pca_mat, centers = k, nstart = 25)
              sil <- silhouette(km$cluster, dist(pca_mat))
              mean(sil[, 3])
            })
            
            
            row_result <- c(
              rank_transformation = rank_trans,
              trans_type = trans_type,
              dist_cal_type = dist_type,
              ord_calc_method = ord_method,
              spearman_rho = rho,
              var_expl_axis1 = var1,
              var_expl_axis2 = var2,
              nmds_stress = stress
            )
            sil_cols <- set_names(as.list(sil_scores), paste0("sil_k", k_range))
            results_list[[counter]] <- c(row_result, sil_cols)
            counter <- counter + 1
            
            
          } )#, silent = TRUE)
        }
      }
    }
  }
  
  # Get all possible columns
  all_columns <- c("rank_transformation", "trans_type", "dist_cal_type", "ord_calc_method",
                   "spearman_rho", "var_expl_axis1", "var_expl_axis2", "nmds_stress",
                   paste0("sil_k", k_range))
  
  # Fill missing columns with NA for incomplete rows
  results_df <- results_list %>%
    map(~{
      missing <- setdiff(all_columns, names(.x))
      if (length(missing) > 0) .x[missing] <- NA
      .x[all_columns]
    }) %>%
    bind_rows() %>%
    mutate(across(all_of(setdiff(all_columns, c("rank_transformation", "trans_type", "dist_cal_type", "ord_calc_method"))), as.numeric))
  
  saveRDS(results_df, result_file)
  return(results_df)
}





#Plot silhouette based on k-value
plot_kmeans_silhouette_scores <- function(scores_df,
                                          component_num = 3,
                                          k_range = 2:10) {
  library(ggplot2)
  library(cluster)
  library(dplyr)
  
  # Extract PCA axes
  pca_matrix <- scores_df[, paste0("Axis", 1:component_num)]
  
  silhouette_scores <- data.frame(
    k = integer(),
    avg_silhouette = numeric()
  )
  
  for (k in k_range) {
    set.seed(42)
    kmeans_res <- kmeans(pca_matrix, centers = k, nstart = 25)
    
    # Silhouette score
    sil <- silhouette(kmeans_res$cluster, dist(pca_matrix))
    avg_sil <- mean(sil[, 3])
    
    silhouette_scores <- rbind(silhouette_scores, data.frame(k = k, avg_silhouette = avg_sil))
  }
  
  # Plot silhouette scores
  sil_plot <- ggplot(silhouette_scores, aes(x = k, y = avg_silhouette)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Average Silhouette Score vs Number of Clusters",
         x = "Number of Clusters (k)",
         y = "Average Silhouette Score") +
    theme_bw()
  
  return(sil_plot)
}


#K-means clustering on PCA or PCoA
analyze_kmeans_on_pca <- function(scores_df,
                                  k_clusters = 2,
                                  component_num = 3,
                                  colors_list = NULL) {
  library(ggplot2)
  library(dplyr)
  library(cluster)
  
  # Extract PCA axes
  pca_matrix <- scores_df[, paste0("Axis", 1:component_num)]
  
  # K-means clustering
  set.seed(42)
  kmeans_res <- kmeans(pca_matrix, centers = k_clusters, nstart = 25)
  scores_df$Cluster <- factor(kmeans_res$cluster)
  
  # PCA plot colored by Type and shaped by Cluster
  pca_plot <- ggplot(scores_df, aes(x = Axis1, y = Axis2, color = Type, shape = Cluster)) +
    geom_point(size = 2) +
    labs(title = paste("K-means Clustering (k =", k_clusters, ") on PCA"),
         x = "PC1", y = "PC2") +
    theme_minimal()
  
  if (!is.null(colors_list)) {
    pca_plot <- pca_plot +
      scale_fill_manual(values = colors_list) +
      scale_color_manual(values = colors_list)
  }
  
  # Bar plot: number of samples per cluster and type
  cluster_counts <- scores_df %>%
    dplyr::count(Cluster, Type)
  
  bar_plot <- ggplot(cluster_counts, aes(x = Cluster, y = n, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("Sample Counts per Cluster by Type (k =", k_clusters, ")"),
         x = "Cluster",
         y = "Number of Samples") +
    theme_minimal()
  
  if (!is.null(colors_list)) {
    bar_plot <- bar_plot +
      scale_fill_manual(values = colors_list)
  }
  
  # Silhouette score
  sil <- silhouette(kmeans_res$cluster, dist(pca_matrix))
  avg_silhouette <- mean(sil[, 3])
  
  return(list(
    scores = scores_df,
    kmeans = kmeans_res,
    silhouette_score = avg_silhouette,
    pca_plot = pca_plot,
    bar_plot = bar_plot
  ))
}


# LDA Data prep functions

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

