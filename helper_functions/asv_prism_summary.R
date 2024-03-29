library(readr)
library(dplyr)

# Check if there are enough command-line arguments
if (length(commandArgs(trailingOnly = TRUE)) < 4) {
  stop("Usage: Rscript asv_prism_summary.R <file_path> <map_file> <rating_column> <treatments> <number_of_taxa_desired>")
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract the column and treatments
file_path <- args[1]
map_file <- args[2]
column <- args[3]
treatments <- unlist(strsplit(args[4], ","))
num <- args[5]

# read in data file
df <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "~", skip = 1)
colnames(df)[1] <- "ASV_ID"
# Check if the column "taxonomy" exists in the data frame
if ("taxonomy" %in% colnames(df)) {
  # Concatenate "taxonomy" with "ASV_ID"
  df$taxonomy <- paste(df$taxonomy, df$ASV_ID, sep = "_")
  df$ASV_ID <- NULL
} else {
  # Rename "ASV_ID" to "taxonomy"
  names(df)[names(df) == "ASV_ID"] <- "taxonomy"
}

# Read map file using readr functions
map <- read_tsv(map_file, comment = "~") %>%
  rename(sampleID = `#SampleID`) %>%
  select({{ column }}, sampleID)

merged_data_list <- list()
for (treatment in treatments) {
  # Filter map dataframe for the current treatment
  sample_columns <- map[map[[1]] == treatment, 2][[1]]
  
  # Filter df dataframe for sample columns corresponding to the current treatment
  df_subset <- df[, c("taxonomy", sample_columns), drop = FALSE]
  
  # Add all sample columns together into one column
  df_subset$summed_values <- rowSums(df_subset[, sample_columns, drop = FALSE], na.rm = TRUE)
  
  # Remove the individual sample columns
  df_subset2 <- df_subset[, !(names(df_subset) %in% sample_columns)]

  colnames(df_subset2)[2] <- treatment 
  
  # Store the aggregated dataframe in the merged_data_list
  merged_data_list[[treatment]] <- df_subset2
}

# Combine all dataframes in merged_data_list by cbind
combined_df <- Reduce(function(x, y) cbind(x, y), merged_data_list)

# Remove duplicate taxonomy column
combined_df <- combined_df[, !duplicated(names(combined_df))]

# List of ambiguous taxa
patterns <- c(
  "k__Bacteria;Other",
  "k__Fungi;Other",
  "k__Eukaryota;Other",
  "k__Bacteria;p__unclassified_Bacteria",
  "k__Fungi;p__unclassified_Fungi",
  "k__Eukaryota;p__unclassified_Eukaryota",
  "k__Bacteria_OR_k__unclassified_;Other",
  "k__Fungi_OR_k__unclassified_;Other",
  "k__Eukaryota_OR_k__unclassified_;Other",
  "k__Unassigned;Other"
)
combined_pattern <- paste(patterns, collapse = "|")

top_taxa <- vector("list", length = num)
counter <- 0
row_counter <- 1

while (counter < num) {  
  # Iterate over each column
  for (col in names(combined_df)[-1]) {  # Exclude the first column which is the taxonomy column
    # Extract column
    counts <- combined_df[[col]]
    
    sorted_indices <- order(counts, decreasing = TRUE) # Sort column in descending order
    
    # Get the corresponding taxa for the current row
    taxa <- combined_df$taxonomy[row_counter]
    
    # Check if the taxa is already in the top_taxa list
    if (!(taxa %in% top_taxa) && !grepl(combined_pattern, taxa)) {
      # If not and taxa does not match the partial pattern, add it to the list
      counter <- counter + 1
      top_taxa[[counter]] <- taxa
    }
    
    # Check if we have found the top num taxa
    if (counter >= num) {
      break  # Exit the loop if we have found the top num taxa
    }
  }
  # Increment row_counter to move to the next row
  row_counter <- row_counter + 1
}

# Filter rows matching top_taxa
top_taxa_df <- combined_df[combined_df$taxonomy %in% top_taxa, ]

# Identify rows not matching top_taxa
other_rows <- combined_df[!(combined_df$taxonomy %in% top_taxa), ]

# Create a new row for "Other" by summing all other rows
other_row <- c("Other", colSums(other_rows[, -1]))  # Assuming numerical columns, adjust if needed

# Combine into a new dataframe
new_df <- rbind(top_taxa_df, other_row)

# Set column names
colnames(new_df) <- colnames(combined_df)

# Clean up taxonomy column
remove_last_semicolon <- function(string) {
  sub("^.+;\\s*", "", string)
}

# Apply the function to the "taxonomy" column of the new dataframe
new_df$taxonomy <- sapply(new_df$taxonomy, remove_last_semicolon)

final <- t(new_df)
colnames(final) <- final[1,]
final <- final[-1,]

write.table(final, file = "prism_input_file.csv", sep = "\t", quote = FALSE, row.names = FALSE)