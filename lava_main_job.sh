#!/bin/bash
# This script processes a SINGLE phenotype pair, parallelizing over loci chunks.
# It's called by the outer GNU Parallel.

module load 2022 
module load parallel/20220722-GCCcore-11.3.0 

# --- Arguments passed from outer GNU Parallel ---
p1="${1}"
p2="${2}"
NUM_INNER_PARALLEL_JOBS="${3:-4}" # Default to 4 if $3 is empty (good practice)

if [ -z "${p1}" ] || [ -z "${p2}" ]; then
    echo "Error [lava_main_job.sh]: Phenotype 1 ('${p1}') and/or Phenotype 2 ('${p2}') not provided or empty."
    exit 1
fi

echo "--- [Pair: ${p1}.${p2}] Starting processing ---"
echo "--- [Pair: ${p1}.${p2}] Host: $(hostname) ---"
echo "--- [Pair: ${p1}.${p2}] Received NUM_INNER_PARALLEL_JOBS = $NUM_INNER_PARALLEL_JOBS ---"

# --- Paths ---
base_dir="/path/to/lava/folder"
results_base_dir="${base_dir}/lava_results"
loci_list_1_100_parallel="/path/to/ref/loci/loci100files/loci_list_1_100_parallel"
lava_script="${base_dir}/lava_scripts/lava_script.sh" 

# Define the specific output directory for this pair's results and logs
output_dir="${results_base_dir}/${p1}/${p1}.${p2}"
mkdir -p "$output_dir" # Ensure it exists, otherwise it will be created

echo "--- [Pair: ${p1}.${p2}] Inner GNU Parallel will run $NUM_INNER_PARALLEL_JOBS loci chunks concurrently. ---"
echo "--- [Pair: ${p1}.${p2}] Output and logs will be in: $output_dir ---"

# --- Check if final aggregated output for this pair already exists ---
final_bivar_output_for_pair_check="${output_dir}/${p1}.${p2}.combined.bivar.lava"
if [[ -f "$final_bivar_output_for_pair_check" ]]; then
    echo "--- [Pair: ${p1}.${p2}] Aggregated bivariate output file '$final_bivar_output_for_pair_check' already exists. Skipping. ---"
    exit 0 # Successful exit, tells outer parallel this job is done.
fi

# Change to the pair's output directory for inner parallel logs and relative paths for results
cd "$output_dir"
if [ $? -ne 0 ]; then
    echo "--- [Pair: ${p1}.${p2}] CRITICAL ERROR: Failed to cd to $output_dir. Exiting. ---"
    exit 1
fi
echo "--- [Pair: ${p1}.${p2}] Current working directory set to: $(pwd) ---"

# --- INNER GNU Parallel for loci chunks ---
echo "--- [Pair: ${p1}.${p2}] Starting INNER GNU Parallel for loci chunks... ---"
cat "$loci_list_1_100_parallel" | \
    parallel \
        --jobs "$NUM_INNER_PARALLEL_JOBS" \
        --eta \
        --joblog "lava_inner.${p1}.${p2}.log" \
        sh "$lava_script" {} "${p1}" "${p2}"

echo "--- [Pair: ${p1}.${p2}] INNER GNU Parallel finished for loci chunks. ---"

# --- Aggregate results (runs in $output_dir, so filenames are relative) ---
echo "--- [Pair: ${p1}.${p2}] Aggregating results... ---"
# final_univ_output="${p1}.${p2}.combined.univ.lava" # Relative filename
final_bivar_output="${p1}.${p2}.combined.bivar.lava" # Relative filename

# # Pattern for chunked files created by lava_script.sh: P1.P2.loci_subset_filename
# chunk_pattern_univ="${p1}.${p2}.*.univ.lava" 
# first_univ_file=$(ls -1 $chunk_pattern_univ 2>/dev/null | head -n 1)
# if [ -n "$first_univ_file" ]; then
#     head -n 1 "$first_univ_file" > "$final_univ_output"
#     find . -maxdepth 1 -type f -name "$chunk_pattern_univ" -print0 | xargs -0 -I {} tail -n +2 "{}" >> "$final_univ_output"
#         # Use find for robustness, especially if many files. Ensure find searches only in current directory for these chunk files.
#     echo "--- [Pair: ${p1}.${p2}] Aggregated univariate results to: $(pwd)/$final_univ_output ---"
# else
#     echo "--- [Pair: ${p1}.${p2}] Warning: No .univ.lava files found matching pattern '$chunk_pattern_univ' in $(pwd) to aggregate. ---"
# fi

chunk_pattern_bivar="${p1}.${p2}.locfile.*.bivar.lava" 
first_bivar_file=$(ls -1 $chunk_pattern_bivar 2>/dev/null | head -n 1)
if [ -n "$first_bivar_file" ]; then
    head -n 1 "$first_bivar_file" > "$final_bivar_output"
    find . -maxdepth 1 -type f -name "$chunk_pattern_bivar" -print0 | xargs -0 -I {} tail -n +2 "{}" >> "$final_bivar_output"
    echo "--- [Pair: ${p1}.${p2}] Aggregated bivariate results to: $(pwd)/$final_bivar_output ---"
else
    echo "--- [Pair: ${p1}.${p2}] Warning: No .bivar.lava files found matching pattern '$chunk_pattern_bivar' in $(pwd) to aggregate. ---"
fi

echo "--- [Pair: ${p1}.${p2}] LAVA analysis completed. ---"
