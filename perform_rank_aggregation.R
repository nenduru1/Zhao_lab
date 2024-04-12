#' Perform Rank Aggregation
#'
#' This function performs rank aggregation using various methods and generates an aggregated ranks CSV file.
#'
#' @param data_path The path to the data file.
#' @param file_prefix The prefix for the output file name.
#' @param ranks_column The column containing ranks in the data.
#' @return A dataframe containing the aggregated ranks.
#' @export


perform_rank_aggregation <- function(data_path, file_prefix, ranks_column) {
  
  # Read data from the specified file
  data <- read.csv(data_path)  # Adjust read function based on your file type
  
  # Compute mean reciprocal rank (MRR)
  # Compute inverse ranks
  data_InverseRank <- data %>%
    mutate(
      inv_ranks = 1 / {{ranks_column}}
    )
  # Reshape data (weights and ID) to wide format (missing weights will be 0)
  data_wide <- data_InverseRank %>%
    dplyr::select(ID, SNP_coord, inv_ranks) %>%
    tidyr::spread(key = ID, value = inv_ranks, fill = 0)
  
  # Calculate average of all columns for each row and add a new column "avg_inv_rank"
  data_wide_MRR <- data_wide %>%
    dplyr::rowwise() %>%
    dplyr::mutate(avg_inv_rank = mean(c_across(-SNP_coord)))
  # Subset MRR and SNP_coord
  data_MRR <- data_wide_MRR %>%
    dplyr::select(SNP_coord, avg_inv_rank)
  
  # import pgs metadata from file. Download excel file from: https://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/
  
  pgs_all_metadata <- read_excel("/data2/gushijimamwesigwa/projects/chat_GPT/pgs_all_metadata.xlsx", sheet = "Score Development Samples")
  # Find the column index by name
  column_index <- which(names(pgs_all_metadata) == "Polygenic Score (PGS) ID")
  # Rename the column
  names(pgs_all_metadata)[column_index] <- "ID"
  
  # Grouping and summarizing the pgs_all_metadata dataframe to obtain the sum for duplicate IDs
  pgs_sample_size <- pgs_all_metadata %>%
    group_by(ID) %>%
    summarize(Sum_Number_of_Individuals = sum(`Number of Individuals`, na.rm = TRUE))
  
  
  # Get unique values of column "ID" from 'data'
  unique_ids <- unique(data$ID)
  
  # Create a new data frame 'sample_size_weight' by left joining 'devpt_sample_size' with unique 'ID'
  sample_size_weight <- data.frame(ID = unique_ids)
  sample_size_weight <- merge(sample_size_weight, pgs_sample_size, by = "ID", all.x = TRUE)
  # compute convex weights
  sample_size_weight$convex_weights <- sample_size_weight$Sum_Number_of_Individuals / sum(sample_size_weight$Sum_Number_of_Individuals)
  
  # Extract weights as a vector
  weights <- sample_size_weight$convex_weights
  weights <- as.numeric(weights)
  
  # Convert data_wide_MRR to a matrix and transpose it
  # Convert data frame to matrix with first column as row names
  data_wide <- as.data.frame(data_wide)
  data_matrix <- data_wide[, -1] # Exclude the first and second column when converting to matrix
  data_matrix <- as.matrix(data_matrix)
  rownames(data_matrix) <- data_wide[, 1]    # Set row names using the values from the first column
  
  # Multiply weights by the transposed data matrix
  weighted_values <- t(t(data_matrix) * weights)
  
  # Compute the weighted mean for each row
  weighted_means <- rowSums(weighted_values)
  weighted_MRR <- data.frame(weighted_means)
  # Add row SNP_coord as a new column
  weighted_MRR$SNP_coord <- row.names(weighted_MRR)
  
  
  #                              *** Using RobustAggreg Package ***
  
  # Reshape data (weights and ID) to wide format (missing weights will be NA)
  data_wide <- data %>% dplyr::select(ID, SNP_coord, {{ranks_column}})
  data_wide <- tidyr::spread(data_wide, key = ID, value = {{ranks_column}}, fill = NA)
  
  # get cut off  values for the RobustRankAggreg
  # Get unique groups in the "ID" column and sort them
  unique_sorted_groups <- sort(unique(data$ID))
  # Count the occurrences of each group in the "ID" column
  group_counts <- table(data$ID)
  # Create a vector of counts in the sorted order of groups
  sorted_counts <- unname(group_counts[match(unique_sorted_groups, names(group_counts))])
  
  # Replace NA weights with maximum rank
  max_rank = nrow(data_wide)
  data_wide <- data_wide %>%
    mutate_all(~ifelse(is.na(.), max_rank, .))
  
  # Convert all columns except the first one to numeric
  data_wide[, -1] <- lapply(data_wide[, -1], as.numeric)
  
  r = as.matrix(data_wide[, -1])
  row.names(r) <- data_wide$SNP_coord
  
  # Apply Min-Max scaling to each column using rescale
  # Load the scales package if not already loaded
  if (!require(scales)) {
    install.packages("scales")
    library(scales)
  }
  
  
  # Apply Min-Max scaling to each column using apply
  r_norm <- as.data.frame(apply(r, 2, rescale))
  r_norm <- as.matrix(r_norm)
  # Print the normalized data
  # print(r_norm)
  
  
  # Use RobustRankAggreg to compute ranks : method options include 'min', 'geom.mean', 'mean', 'median', 'stuart' or 'RRA'
  # aggreg_rank <- aggregateRanks(rmat = r, method = "geom.mean", exact = TRUE, topCutoff = sorted_counts)
  aggreg_rank_RRA <- aggregateRanks(rmat = r_norm, method = "RRA", exact = TRUE)
  aggreg_rank_stuart <- aggregateRanks(rmat = r_norm, method = "stuart", exact = TRUE, topCutoff = sorted_counts)
  aggreg_rank_min <- aggregateRanks(rmat = r_norm, method = "min", exact = TRUE, topCutoff = sorted_counts)
  aggreg_rank_geo <- aggregateRanks(rmat = r_norm, method = "geom.mean", exact = TRUE, topCutoff = sorted_counts)
  aggreg_rank_mean <- aggregateRanks(rmat = r_norm, method = "mean", exact = TRUE, topCutoff = sorted_counts)
  aggreg_rank_median <- aggregateRanks(rmat = r_norm, method = "median", exact = TRUE, topCutoff = sorted_counts)
  
  
  # Function to Replace "Score" and "Name" in the dataframe
  rename_columns <- function(dataframe, substitution) {
    # Rename "Name" to "SNP_coord" and substitute "Score" with the provided value
    colnames(dataframe) <- sub("^Name$", "SNP_coord", colnames(dataframe))
    colnames(dataframe) <- sub("^Score$", substitution, colnames(dataframe))
    
    return(dataframe)
  }
  
  
  # Replace "Score" and "Name" in the dataframe
  aggreg_rank_RRA <- rename_columns(aggreg_rank_RRA, "RRA_score")
  aggreg_rank_stuart <- rename_columns(aggreg_rank_stuart, "Stuart_score")
  aggreg_rank_min <- rename_columns(aggreg_rank_min, "Min_score")
  aggreg_rank_geo <- rename_columns(aggreg_rank_geo, "Geo_score")
  aggreg_rank_mean <- rename_columns(aggreg_rank_mean, "Mean_score")
  aggreg_rank_median <- rename_columns(aggreg_rank_median, "Med_score")
  
  #reshape data to wide dataframe with gene annotation
  # Reshape data (weights and ID) to wide format (missing weights will be 0)
  data_wide_anno <- data %>%
    dplyr::select(ID, SNP_coord, hm_chr, hm_pos, hm_rsID, annot.symbol, annot.type, ranks) %>%
    tidyr::spread(key = ID, value = ranks, fill = NA)
  
  # Add MRR ranks
  data_wide_anno_MRR <- left_join(data_wide_anno, data_MRR, by = "SNP_coord")
  data_wide_anno_MRR <- left_join(data_wide_anno_MRR, weighted_MRR, by = "SNP_coord")
  
  
  # add RobustRankAggreg scores to data_wide_anno 
  data_aggreg_rank <- left_join(data_wide_anno_MRR, aggreg_rank_RRA, by = "SNP_coord") %>%
    left_join(aggreg_rank_stuart, by = "SNP_coord") %>% 
    left_join(aggreg_rank_min, by = "SNP_coord") %>% 
    left_join(aggreg_rank_geo, by = "SNP_coord") %>% 
    left_join(aggreg_rank_mean, by = "SNP_coord") %>% 
    left_join(aggreg_rank_median, by = "SNP_coord")
  
  # Rank columns such that top rank is 1 ...
  data_aggreg_rank <- data_aggreg_rank %>%
    mutate(
      RRA_score = rank(RRA_score),
      Stuart_score = rank(Stuart_score),
      Min_score = rank(Min_score),
      Geo_score = rank(Geo_score),
      Mean_score = rank(Mean_score),
      Med_score = rank(Med_score),
      avg_inv_rank = rank(-avg_inv_rank),
      weighted_means = rank(-weighted_means)
    )
  
  # Final dataframe with aggregated ranks
  final_data <- data_aggreg_rank %>% dplyr::select(hm_rsID, hm_chr, hm_pos, annot.symbol, annot.type, weighted_means, avg_inv_rank, RRA_score, Stuart_score, Min_score, Geo_score, Mean_score, Med_score)
  
  # Save the final dataframe as a CSV file
  file_name <- paste0(file_prefix, "_aggregated_ranks.csv")
  write.csv(final_data, file = file_name, row.names = FALSE)
  
  # Return the final aggregated dataframe
  return(final_data)
}
