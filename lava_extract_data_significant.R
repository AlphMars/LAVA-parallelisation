# Script for Bonferroni/FDR significant loci extraction
# Called by the command in 'lava_extract_data.sh'

#################################################
# ===< 0. EXTRACT ARGUMENT AND PREPARE ENV >=== #
#################################################
arg = commandArgs(T)
wd = arg[1]
# files = arg[2] # in case you want to provide a list as argument instead of using the method below

if (!"dplyr" %in% installed.packages()) {
  install.packages("dplyr")}
if (!"fs" %in% installed.packages()) {
  install.packages("fs")}
library(dplyr)
library(fs)
setwd(wd)

####################################
# ===< 1. PREPARE INPUT FILES >=== #
####################################
files <- fs::dir_ls(path = ".", glob = "*/*/*results.bivar.lava", recurse = TRUE)
  # find files matching the pattern (recursive search)
files <- as.character(files)
  # convert file paths to character strings
tables <- lapply(files, read.table, header = TRUE)
  # read in all files as separate tables, store in a list
table_combined <- bind_rows(tables, .id="file") 
  # combine all tables into a single table, removing duplicate headers
    # '.id = "file"' : a new column called "file" is created that indicates which original table each row came from.
table_combined_unique <- distinct(table_combined, .keep_all = TRUE) 
  # Remove duplicate rows if any and keep the first occurrence
    # '.keep_all=T' : makes sure that all columns are preserved in the final output

#############################################
# ===< 2. CALCULATE BONFERRONI'S ALPHA >=== #
#############################################
num_rows <- nrow(table_combined_unique)
bonf_alpha <- 0.05/num_rows

table_corrected_bonf <- table_combined_unique %>% filter(p < bonf_alpha)
  # filter bonf significant regions
table_nominal_sign <- table_combined_unique %>% filter(p < 0.05)
  # filter nominal significant regions

# write combined table to a text file with headers
write.table(table_combined_unique, file = "./all_pheno_pairs.all_loci.results.bivar.lava", sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(table_corrected_bonf, file = paste0("./all_pheno_pairs.all_loci.results.bivar.lava",".bonf_sign_",bonf_alpha), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(table_nominal_sign, file = paste0("./all_pheno_pairs.all_loci.results.bivar.lava",".nominal_sign"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

################################################
# ===< 3. APPLY FDR CORRECTION FOR q=0.05 >=== #
################################################
table_combined_unique$q_value <- p.adjust(table_combined_unique$p, method = "fdr", n = length(table_combined_unique$p))
  # calculate q value
table_corrected_fdr <- table_combined_unique[table_combined_unique$q_value <= 0.05, ]
  # filter for q-value <= 0.05
# write combined table to a text file with headers
write.table(table_corrected_fdr, file = paste0("./all_pheno_pairs.all_loci.results.bivar.lava",".fdr_sign"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
# Split the data by the "file" column
data_list_fdr_sign_and_not <- split(table_combined_unique, table_combined_unique$file) # also not FDR significant
data_list_fdr <- split(table_corrected_fdr, table_corrected_fdr$file)
data_list_bonf <- split(table_corrected_bonf, table_corrected_bonf$file)

################################
# ===< 4. FILTERED OUTPUT >=== #
################################
# ===< A. Loop through each df in the list and write to a separate file all regions with q-values >=== #
for (i in 1:length(data_list_fdr_sign_and_not)) {
  # Extract phenotypes names in the dataframe
  phen1_name <- unique(data_list_fdr_sign_and_not[[i]]$phen1)
  phen2_name <- unique(data_list_fdr_sign_and_not[[i]]$phen2)
  # Indicate the output directory and set the flag to determine whether it exists
  output_dir1 <- paste0("./", phen1_name, "/", phen1_name, ".", phen2_name)
  output_dir2 <- paste0("./", phen2_name, "/", phen2_name, ".", phen1_name)
  dir_exists <- file.exists(output_dir1) || file.exists(output_dir2)
  # Choose the output directory based on which one exists
  if (file.exists(output_dir1)) {
    output_dir <- output_dir1
  } else if (file.exists(output_dir2)) {
    output_dir <- output_dir2
  } else {
    # If neither directory exists, print a message and move on to the next iteration
    message(paste("Neither", output_dir1, "nor", output_dir2, "exist. Skipping", phen1_name, "-", phen2_name))
    next
  }
  # Create the output filename
  output_filename <- paste0(phen1_name, ".", phen2_name, ".all_loci.results.bivar.lava.fdr_qvalues")
  # Write the dataframe to a file with the new filename
  write.table(data_list_fdr_sign_and_not[[i]], file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ===< B. Loop through each df in the list and write to a separate file the FDR-significant regions >=== #
for (i in 1:length(data_list_fdr)) {
  # Extract phenotypes names in the dataframe
  phen1_name <- unique(data_list_fdr[[i]]$phen1)
  phen2_name <- unique(data_list_fdr[[i]]$phen2)
  # Indicate the output directory and set the flag to determine whether it exists
  output_dir1 <- paste0("./", phen1_name, "/", phen1_name, ".", phen2_name)
  output_dir2 <- paste0("./", phen2_name, "/", phen2_name, ".", phen1_name)
  dir_exists <- file.exists(output_dir1) || file.exists(output_dir2)
  # Choose the output directory based on which one exists
  if (file.exists(output_dir1)) {
    output_dir <- output_dir1
  } else if (file.exists(output_dir2)) {
    output_dir <- output_dir2
  } else {
    # If neither directory exists, print a message and move on to the next iteration
    message(paste("Neither", output_dir1, "nor", output_dir2, "exist. Skipping", phen1_name, "-", phen2_name))
    next
  }
  # Create the output filename
  output_filename <- paste0(phen1_name, ".", phen2_name, ".all_loci.results.bivar.lava.fdr_sign")
  # Write the dataframe to a file with the new filename
  write.table(data_list_fdr[[i]], file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  # Select the desired columns from the DataFrame and write to a file with the new filename
  selected_data <- data_list_fdr[[i]] %>%
  dplyr::select(chr, start, stop)
  write.table(selected_data, file = paste0(output_filename,"_for_biomart"), sep = "\t", quote = FALSE, row.names = FALSE)
  # with cromosomal regions according biomart format
  selected_data_b <- data_list_fdr[[i]] %>%
  dplyr::select(chr, start, stop) %>%
  tidyr::unite("chromosomal_region", chr, start, stop, sep=":")
  write.table(selected_data_b, file = paste0(output_filename,"_for_biomart_chromosomal_region"), sep = "\t", quote = FALSE, row.names = FALSE)
  # create bed files (e.g. for loci2path) starting from LAVA correlated regions for each couple of phenos (tab delimited chr1 startpb endpb with no headers)
  bed <- data_list_fdr[[i]] %>% 
  dplyr::select(chr, start, stop) %>% 
  mutate(chr = paste0("chr", chr))
  dir.create("./lava_bed", showWarnings = FALSE, recursive = TRUE) # The recursive = TRUE argument ensures that any missing parent directories in the specified path are also created.
  write.table(bed, file = paste0("./lava_bed/", phen1_name, ".", phen2_name, ".all_loci.results.bivar.lava.fdr_sign.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ===< C. Loop through each dataframe in the list and write to a separate file the BONFERRONI significant regions >=== #
# NB: some pairs may not have bonf significant regions!!!
for (i in 1:length(data_list_bonf)) {
  # Extract phenotypes names in the dataframe
  phen1_name <- unique(data_list_bonf[[i]]$phen1)
  phen2_name <- unique(data_list_bonf[[i]]$phen2)
  # Indicate the output directory and set the flag to determine whether it exists
  output_dir1 <- paste0("./", phen1_name, "/", phen1_name, ".", phen2_name)
  output_dir2 <- paste0("./", phen2_name, "/", phen2_name, ".", phen1_name)
  dir_exists <- file.exists(output_dir1) || file.exists(output_dir2)
  # Choose the output directory based on which one exists
  if (file.exists(output_dir1)) {
    output_dir <- output_dir1
  } else if (file.exists(output_dir2)) {
    output_dir <- output_dir2
  } else {
    # If neither directory exists, print a message and move on to the next iteration
    message(paste("Neither", output_dir1, "nor", output_dir2, "exist. Skipping", phen1_name, "-", phen2_name))
    next
  }
  # Create the output filename
  output_filename <- paste0(phen1_name, ".", phen2_name, ".all_loci.results.bivar.lava.bonf_sign_",bonf_alpha)
  # Write the dataframe to a file with the new filename
  write.table(data_list_bonf[[i]], file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  # Select the desired columns from the DataFrame and write to a file with the new filename
  selected_data <- data_list_bonf[[i]] %>% dplyr::select(chr, start, stop)
  write.table(selected_data, file = paste0(output_filename,"_for_biomart"), sep = "\t", quote = FALSE, row.names = FALSE)
  # create bed files (e.g. for loci2path) starting from LAVA correlated regions for each couple of phenos (tab delimited chr1 startpb endpb with no headers)
  bed_bonf <- data_list_bonf[[i]] %>% 
  dplyr::select(chr, start, stop) %>% 
  mutate(chr = paste0("chr", chr))
  dir.create("./lava_bed", showWarnings = FALSE, recursive = TRUE) # The recursive = TRUE argument ensures that any missing parent directories in the specified path are also created.
  write.table(bed, file = paste0("./lava_bed/", phen1_name, ".", phen2_name, ".all_loci.results.bivar.lava.bonf_sign_",bonf_alpha, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
