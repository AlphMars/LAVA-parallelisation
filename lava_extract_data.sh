#!/bin/bash
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --job-name=lava_results
#SBATCH --output=/path/to/lava/folder/lava_scripts/lava.%A.out 
#SBATCH --error=/path/to/lava/folder/lava_scripts/lava.%A.err 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=...

module load 2022
module load R/4.2.1-foss-2022a


base_dir="/path/to/lava/folder"
pheno_pairs_file="${base_dir}/lava_pheno_pairs.tsv"
results_base_dir="${base_dir}/lava_results"
n_pairs=$(wc -l < "${pheno_pairs_file}") # Use < for wc -l (faster)

#########################
######## Tables #########
#########################
for i in $(seq 1 "${n_pairs}"); do
    p1=$(awk 'NR=='$i' {print $1}' "${pheno_pairs_file}")
    p2=$(awk 'NR=='$i' {print $2}' "${pheno_pairs_file}")

    cd ${results_base_dir}/${p1}/${p1}.${p2}/
        # Navigate to the directory containing LAVA outputs
    grep -E "${p1} ${p2} | ${p2} ${p1}" *bivar* | sort -gk15 > "${p1}.${p2}.all_loci.results.bivar.lava"
        # Search through all the *bivar* files for lines matching the pattern:
            # pheno1 pheno2 OR pheno2 pheno1
            # '|' is the OR operator in grep
            # -E is the extended regular expression option, allowing for more complex patterns. If you don't need it, you'd want to escape the '|' operator like this: '\|'
        # sort the output by number 
            # '-g', allows for general numeric format, including scientific notation
            # '-k15' sort based on the 15th column 
    sed -i -e 's/.*://g' ${p1}.${p2}.all_loci.results.bivar.lava 
        # '-e' for the regex, 
            # '.' matches any character
            # '*' repeats this any number of times (including zero)
            # ':' matches a colon
        # it probably removes the `filename` prefix that `grep` adds to each line in the aggregated file
    # echo -e "locus chr start stop n.snps n.pcs phen1 phen2 rho rho.lower rho.upper r2 r2.lower r2.upper p" |\
    # cat - "${p1}.${p2}.all_loci.results.bivar.lava" > "${p1}.${p2}.all_loci.results.bivar.lava.def"
    headers="locus chr start stop n.snps n.pcs phen1 phen2 rho rho.lower rho.upper r2 r2.lower r2.upper p"
    final_file="${p1}.${p2}.all_loci.results.bivar.lava"
    sed -i "1i $headers" "$final_file"
        # '1' : addresses the first line of the file.
        # 'i' : insert command, it inserts the text that follows before the addressed line.
    echo "Done with ${p1}.${p2}"
done


#################################################
#### Extract Bonferroni/FDR significant loci ####
#################################################

# # If you want, you can specify only a subset of files on which to run the R script:
# files=( "" )

# # uncomment \ and ${files} at the end of the Rscript line below to use this option
# # remember to do the same in the R script!

Rscript lava_extrat_data_significant.R \
    ${results_base_dir} # \
    # ${files}
