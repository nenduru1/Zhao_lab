#' Generate full dataset
#'
#' This function merges all PGS score files into one dataset.
#'
#' @param trait_term The trait term to retrieve PGS scores for.
#' @param directory_path The directory path where the PGS score files are located.
#' @param output_file_prefix The prefix for the output file name.
#' @return The full dataset.
#' @export


# ****        generate_full_dataset function to merge all PGS score files into one dataset
generate_full_dataset <- function(trait_term, directory_path, output_file_prefix) {
  print(paste("Trait term received:", trait_term))
  library(dplyr)
  library(tidyr)
  library(quincunx)
  # Retrieve PGS traits related to the provided trait term
  PGS_traits <- quincunx::get_traits(trait_term = trait_term, exact_term = FALSE)
  trait_vector <- PGS_traits@pgs_ids[4]
  
  # Define the file suffix
  suffix <- "_hmPOS_GRCh38.txt"
  
  # Read PGS files related to the provided trait from the specified directory
  results <- lapply(trait_vector, function(item) {
    full_path <- file.path(directory_path, paste0(item, suffix))
    read_scoring_file(full_path)
  })
  
  # Extract all file names
  file_names <- names(results$pgs_id)
  
  # Process data and store in a list
  data_list <- lapply(file_names, function(file_name) {
    cleaned_name <- basename(file_name)
    cleaned_name <- sub("_hmPOS_GRCh38.txt$", "", cleaned_name)
    
    data <- results$pgs_id[[file_name]]$data
    data$ID <- cleaned_name
    
    return(data)
  })
  
  # Merge datasets into a single dataframe
  full_dataset <- dplyr::bind_rows(data_list)
  
  # Save the full dataset as a CSV file
  output_file <- paste0(output_file_prefix, "_full_dataset.csv")
  write.csv(full_dataset, file = output_file, row.names = FALSE)
  
  cat("Full dataset saved as", output_file, "\n")
}
