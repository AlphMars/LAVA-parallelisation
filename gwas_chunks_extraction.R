# For each GWAS summary statistic file, extract genomic chunks based on coordinates defined in a LAVA bivariate analyses results file.
# Maybe it would also be possible to do it on the ones from the univariate analyses, even though it would be less consequential

# ==> 1. Load packages <==
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))


# ==> 2. Define configuration and paths <==
chunks_file <- "/path/to/lava/folder/lava_results/all_pheno_pairs.all_loci.results.bivar.lava.fdr_sign.uniq"
output_dir <- "/path/to/fine/mapping/gwas_chunks"

# Create the output directory if it doesn't exist.
# `recursive = TRUE` creates parent directories if needed.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==> 3. Define Input Files <==
# Be careful that the order of the traits corresponds to the order of the sumstats files in the sumstats variable
traits <- c(...)    # e.g. c("trait1", "trait2", "trait3", ...)

sumstats <- c()     # e.g. c("/path/to/gwas/trait1.sumstats", "/path/to/gwas/trait2.sumstats", ...)


# ==> 4. Read Master Chunk File <==
# Read the file once outside the loop/parallel process for efficiency
# Use fread for faster reading
all_chunks <- fread(chunks_file, header = TRUE)

# ==> 5. Define the Processing Function <==
# It will process one trait at a time and looped across all traits in parallel
process_trait <- function(i) {
    # Get the trait name and corresponding sumstat
    current_trait <- traits[i]
    gwas_file <- sumstats[i]

    # Read the full GWAS summary statistics
    # `tryCatch` will prevent the entire script from failing if one file is corrupt or missing
    gwas_data <- tryCatch({
        fread(gwas_file) %>%
        rename_all(tolower) # Standardize column names (lowercase headers)
    }, error = function(e) {
        cat(paste0("ERROR: Could not read or process file for trait '", current_trait, "': ", gwas_file, "\n", e$message, "\n"))
        return(NULL) # Return NULL on error
    })

    # If reading the file failed, stop processing for this trait
    if (is.null(gwas_data)) return(NULL)

    # Armonize column names. Modify as needed.
    # This block makes the script robust to different column naming conventions.
    if ("chromosome" %in% names(gwas_data)) setnames(gwas_data, "chromosome", "chr")
    if ("position" %in% names(gwas_data)) setnames(gwas_data, "position", "bp")
    if ("pos" %in% names(gwas_data)) setnames(gwas_data, "pos", "bp")

    # Filter the master chunks data frame to get loci for the current trait
    trait_chunks <- all_chunks[phenotype == current_trait]

    # Check if there are any chunks to process for this trait
    if (nrow(trait_chunks) == 0) {
        cat(paste0("INFO: No significant chunks found for trait: ", current_trait, "\n"))
        return(NULL) # Nothing to do, so exit function for this trait
    }

    # Loop through each chunk/locus for the current trait
    for (j in 1:nrow(trait_chunks)) {
        locus_info <- trait_chunks[j, ]

    # Filter the GWAS data for SNPs within the current chunk's coordinates
    # Using `locus_info$` is cleaner and avoids variable name confusion
    filtered_data <- gwas_data %>%
      filter(chr == locus_info$chr & bp >= locus_info$start & bp <= locus_info$end)

    # If any SNPs are found in the chunk, write them to a file
    if (nrow(filtered_data) > 0) {
      # Construct a clear output filename
      output_filename <- file.path(output_dir, paste0(current_trait, ".chunk_", locus_info$locus))

      # Write the filtered data to the output file
      fwrite(filtered_data,
             file = output_filename,
             sep = " ",
             quote = FALSE,
             row.names = FALSE,
             col.names = TRUE)
    }
  }
  
  cat(paste0("SUCCESS: Finished processing all chunks for trait: ", current_trait, "\n"))
  # Return the number of chunks processed for this trait for a summary
  return(nrow(trait_chunks))
}


# ==> 6. Run the Processing in Parallel <==
# Determine the number of cores to use.
# `Sys.getenv("SLURM_CPUS_PER_TASK")` reads the value from your sbatch script.
# Fallback to a default of 4 if not running in a Slurm job.
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 4))

# Use mclapply to apply the `process_trait` function to each index from 1 to length(traits).
# It will distribute the work across `num_cores`.
# The loop runs from 1 to length(traits), and the index is passed as `i` to the function.
results <- mclapply(1:length(traits), process_trait, mc.cores = num_cores)

# ==> 7. Final Summary <==
cat("\n\n---------------------------------\n")
cat("Parallel processing complete.\n")
cat("Summary of chunks written per trait:\n")
# Create a data frame for a nice summary printout
summary_df <- data.frame(
    trait = traits,
    chunks_processed = sapply(results, function(x) ifelse(is.null(x), 0, x))
)
print(summary_df)
cat("---------------------------------\n")
cat("Script finished.\n")


# --------------------------------------------------------------------------- #
# ---------------------------------- Notes ---------------------------------- #
# --------------------------------------------------------------------------- #
# Parallelization: It uses parallel::mclapply to process multiple GWAS files, 
    # making full use of the 16 cores you requested. 
    # This will be dramatically faster than a simple for loop.
# Encapsulation: The main logic is wrapped in a function (process_trait).
    # This makes the code cleaner and is necessary for mclapply.
# Robustness:
    # It uses file.path() to construct paths, which is safer than pasting strings with /.
    # It creates the output directory automatically (dir.create).
    # A tryCatch block is included so that if one file fails to load, the entire job doesn't crash.
# Dynamic Core Count: The script automatically detects the number of cores allocated by Slurm 
    # (Sys.getenv("SLURM_CPUS_PER_TASK")), making it adaptable to different sbatch configurations.
# Fast I/O: It uses data.table::fwrite() for writing files, which is significantly faster.