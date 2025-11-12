#!/bin/sh
# This script is called by the INNER GNU Parallel for a single loci_chunk/p1/p2 combination

module load 2022
module load R/4.2.1-foss-2022a
  # You could actually call the module in the main loop job script, it would be inherited anyway
  # However, I find it more reasonable to have it here so that I don't miss it

loci_subset="$1" 
p1="$2"               
p2="$3"               

# The ref_loci path uses $1 (loci_subset) directly
ref_loci="/path/to/ref/loci/loci100files/${loci_subset}"
echo "*** [Chunk: $loci_subset | Pair: ${p1}.${p2}] Starting R script. ***"

base_dir="/path/to/lava/folder" 
lava_R_script="${base_dir}/lava_scripts/lava_script.R"
prefix_ref_genome="/path/to/ref_genotype/g1000_eur_maf005"
input_info="${base_dir}/lava_input_infos.tsv"
overlap="${base_dir}/lava_sample_overlap.tsv"
pheno="${p1};${p2}"
output_prefix="${base_dir}/lava_results/${p1}/${p1}.${p2}/${p1}.${p2}.${loci_subset}"

echo "    R script path: $lava_R_script"
echo "    Ref Loci (chunk): $ref_loci"
echo "    Phenos string: $pheno"
echo "    Output file(s): $output_prefix"

Rscript "$lava_R_script" \
    "$prefix_ref_genome" \
    "$ref_loci" \
    "$input_info" \
    "$overlap" \
    "$pheno" \
    "$output_prefix"

R_exit_code=$?
if [ $R_exit_code -ne 0 ]; then
    echo "*** [Chunk: $loci_subset | Pair: ${p1}.${p2}] Error: R script failed with exit code $R_exit_code. ***"
    # Consider also writing this error to a specific error file for this chunk/pair if needed for easier debugging.
    exit $R_exit_code # Propagate error; GNU Parallel will log this in its joblog.
fi

echo "*** [Chunk: $loci_subset | Pair: ${p1}.${p2}] Finished R script. ***"
