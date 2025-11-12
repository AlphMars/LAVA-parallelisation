#!/bin/bash
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=1-20:00:00
#SBATCH --job-name=lava
#SBATCH --output=/path/to/lava/folder/lava_scripts/lava.%A.out 
#SBATCH --error=/path/to/lava/folder/lava_scripts/lava.%A.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alfonso.martone@radboudumc.nl

echo "----------------------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Running all LAVA phenotype pairs within a single Slurm allocation."
echo "Total Cores allocated: $SLURM_CPUS_PER_TASK" # This will be 32
echo "----------------------------------------------------"

# --- Modules and Environment ---
module load 2022 # Or specific year for your Snellius environment
module load parallel/20220722-GCCcore-11.3.0 # Or latest parallel

# --- Configuration for Parallelism ---
export OMP_NUM_THREADS=1
export NUM_CONCURRENT_PAIRS_CONFIG=3
# How many PHENOTYPE PAIRS to process concurrently by this outer GNU Parallel.
export NUM_INNER_PARALLEL_JOBS_CONFIG=8
# How many LOCI CHUNKS to process concurrently PER PAIR by the inner GNU Parallel.

# Sanity check for resource allocation (conceptual)
    # Total R scripts running = NUM_CONCURRENT_PAIRS_CONFIG * NUM_INNER_PARALLEL_JOBS_CONFIG
    # Total threads = Total R scripts * OMP_NUM_THREADS
    # Example for this script's SBATCH request (32 cores):
    # If NUM_CONCURRENT_PAIRS_CONFIG = 2 and NUM_INNER_PARALLEL_JOBS_CONFIG = 8,
    # then Total R scripts = 2 * 8 = 16 R scripts.
    # Total threads = 16 R scripts * 2 OMP_NUM_THREADS = 32 threads. This matches SBATCH --cpus-per-task=32.
    # The current settings (1 pair * 8 inner jobs * 2 OMP = 16 threads) uses half the 32 allocated cores.
    # This is fine for testing, but for production you'd want to use settings that fully utilize the 32 cores.

# --- Paths ---
base_dir="/path/to/lava/folder"
pheno_pairs_file="${base_dir}/lava_pheno_pairs.tsv"
lava_main_job="${base_dir}/lava_scripts/lava_main_job.sh" 
results_base_dir="${base_dir}/lava_results"

# --- Pre-create all top-level output directories for pairs to avoid race conditions if GNU parallel starts multiple pair scripts very quickly.
echo "Pre-creating output directories..."
n_pairs=$(wc -l < "${pheno_pairs_file}") # Use < for wc -l (faster)
for i in $(seq 1 "${n_pairs}"); do
    p1=$(awk 'NR=='$i' {print $1}' "${pheno_pairs_file}")
    p2=$(awk 'NR=='$i' {print $2}' "${pheno_pairs_file}")
    output_dir="${results_base_dir}/${p1}/${p1}.${p2}"
    mkdir -p "$output_dir"
done
# echo "Output directory pre-creation finished."


echo "Outer GNU Parallel will run $NUM_CONCURRENT_PAIRS_CONFIG phenotype pairs concurrently."
echo "Each pair will internally parallelize $NUM_INNER_PARALLEL_JOBS_CONFIG loci chunks."

tail -n +21 "$pheno_pairs_file" | \
    parallel \
        --jobs "$NUM_CONCURRENT_PAIRS_CONFIG" \
        --colsep '\s+' \
        --eta \
        --joblog "${output_dir}/lava_outer.log" \
        sh "$lava_main_job" {1} {2} "$NUM_INNER_PARALLEL_JOBS_CONFIG" 
        # {1} is p1, {2} is p2, $NUM_INNER_PARALLEL_JOBS_CONFIG is passed as $3

echo "----------------------------------------------------"
echo "All phenotype pairs processing initiated by outer GNU Parallel."
echo "Job completed."
echo "----------------------------------------------------"
