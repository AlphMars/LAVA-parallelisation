#!/bin/bash
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --job-name=gwas_chunks_extraction
#SBATCH --output=/path/to/fine/mapping/gwas_chunks_extraction.%A.out 
#SBATCH --error=/path/to/fine/mapping/gwas_chunks_extraction.%A.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=XXX
#SBATCH --mem=32G

# Load necessary modules
module load 2022
module load R/4.2.1-foss-2022a

# Execute the R script
Rscript gwas_chunks_extraction.R