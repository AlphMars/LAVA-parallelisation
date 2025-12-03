# Local Analysis of [co]Variant Annotation parallelization
Collection of notes and scripts on how to run LAVA on multiple trait pairs in a SLURM-using HPC cluster. 
This tool allows the user to estimate local genetic correlation based on GWAS summary statistics and predefined genomic chunks. Check the original [paper][1] and the [GitHub page][2] for more.

Running LAVA
-----------------------------------
It's possible to run LAVA independently, on only selected sumstats, with no need for parallelization. However, I created a slurm job that can be useful in case you need to run LAVA on multiple traits in the most efficient way. These scripts were combine Giuseppe Fanelli's scripts and the original [LAVA cluster setup][3], perfected thanks to Margo Raijmakers, but even like this it took around 27h to process 110 trait pairs. 

>:memo: Please **note** that these scripts assume the following folder structure:
>
>- "/path/to/lava/folder/" &rarr; folder containing the "lava_pheno_pairs.tsv" file
>- "/path/to/lava/folder/lava_results" &rarr; output folder
>- "/path/to/lava/folder/lava_scripts" &rarr; folder containing the presents scripts
>- "/path/to/ref/loci/loci100files/loci_list_1_100_parallel" &rarr; folder storing the loci chunks
> - "/path/to/ref_genotype" &rarr; path to folder containing the reference genotype (e.g. "g1000_eur_maf005")
> 
> To avoid confusion, I suggest to keep the same structure and to adapt just paths.

> :memo: These scripts load the needed modules, which may need adaptation. 


### "lava_main_job_loop.sh"
This is the master job, the one you need to submit. It will loop the "lava_main_job.sh" through the different trait pairs.

### "lava_main_job.sh"
This script processes a single phenotype pair (provided by the master job), parallelising over loci chunks the inner loop represented by "lava_script.sh".

### "lava_script.sh"
This is the core job, finally running the "lava_script.R" script for each trait pair defined by "lava_main_job_loop.sh" and each locus defined by "lava_main_job.sh". It provides the R script with all the arguments needed to run the analysis on that locus in that specific trait pair.

### "lava_script.R"
Adapted from Giuseppe's script. It runs the analysis on a loci clump specified in the arguments. command line arguments, specifying input/output file names and phenotype subset.
>It shouldn't need any adaption. However, you may need to install LAVA on your R the first time if the installation command fails. Check the [GitHub][2] for more.

### Needed files
#### Sumamry statistics    
- Requested columns:
  - SNP / ID / SNPID_UKB/ SNPID / MarkerName / RSID / RSID_UKB
  - A1 / ALT: effect allele
  - A2 / REF: reference allele
  - N / NMISS / N_analyzed: number of samples
  - Z / T / STAT / Zscore: if provided, no p-values or coefficients are needed; otherwise, please provide both

The qc steps recommended are the same you run, for example, for GenomicSEM. However, you do need to extract the needed columns.

#### Other input files
- **lava_input_infos.tsv** : File summarising the info of the included traits in the following format:

        phenotype	cases	controls	filename
        B_fluid_int	1	0	/path/to/lava/folder/B_fluid_int_baseline_for_lava.tsv
        B_num_mem	1	0	/path/to/lava/folder/B_num_mem_baseline_for_lava.tsv
        B_rev_mean_RT	1	0	/path/to/lava/folder/B_rev_mean_RT_baseline_for_lava.tsv
        B_rev_pairs_match	1	0	/path/to/lava/folder/B_rev_pairs_match_baseline_for_lava.tsv
        I_Matrix	1	0	/path/to/lava/folder/I_Matrix_imaging_for_lava.tsv


- **lava_pheno_pairs.tsv** : Tab-separated list of all the phenotype pairs

- **lava_sample_overlap.tsv** : Sample overlap must be estimated. It can be easily extracted from LDSC results calculated with GenomicSEM, e.g.:
```r
# load LDSC results
load("/path/to/LDSC/results/LDSCoutput_all_pheno.RData")
# calculate the sample overlap fromt the intercepts
corr <- round(cov2cor(LDSCoutput_all_pheno$I), 5)
colnames(corr) = rownames(corr) <- colnames(LDSCoutput_all_pheno$S)

# save the results
write.table(corr, file = "/path/to/lava/folder/lava_sample_overlap.tsv", sep = "\t", col.names = T, row.names = T, quote = F)
```

- **Reference genotype** files in plink format (bed/bim/fam). Download the one that suits you the best (e.g. the g1000_eur). However, remember you may need to filter for MAF>0.05.
  
- **blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile**. Reference loci file dowloaded from the [cluster setup][4]. You can split the file in clumps in R:
```r
# create the directory in which to save the loci files
dir.create("/path/to/ref_loci/loci100files/")
# load the ref loci file
locfile = read.table("/path/to/ref_loci/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile",header=T)
# split the file based on the number of rows (aka, number of loci) devided by 100
splits = split(locfile$LOC, rep(1:100, each=ceiling(nrow(locfile)/100))[1:nrow(locfile)])
# save different loci files for each clump
for (i in names(splits)) { 
    write.table(subset(locfile, LOC %in% unlist(splits[i])), paste0("/path/to/ref_loci/loci100files/locfile.",i), row.names=F, quote=F)
}
```
- **loci_list_1_100_parallel**: list of the loci files
```sh
ls locfile* > /path/to/ref_loci/loci100files/loci_list_1_100_parallel
```


Results extraction
------------------------------------
### lava_extract_data.sh
Creates the input files from LAVA results and calls "lava_extract_data_significant.R" on these.
It's also possible to select a subset of files.

### lava_extract_data_significant.R
R script to extract FDR/Bonferroni significant results. It also creates bed files that can be used for example with loci2path analyses.

Extract GWAS chunks
------------------------------------
Subset the GWAS summary statistics based on the chunks used in LAVA. These can be used for colocalization and fine-mapping analyses (e.g. mvSuSiE).
### gwas_chunks_extraction.sh
Adapt modules and settings as needed.

### gwas_chunks_extraction.R
Adapt as needed (paths and traits).

### Needed files
- **all_pheno_pairs.all_loci.results.bivar.lava.fdr_sign.uniq**: list of all loci significant for each trait from LAVA bivariate analyses

        locus   chr     start   end     phenotype
        9       1       6136815 7711794 B_fluid_int
        17      1       15755231        16732168        B_fluid_int
        58      1       67761891        68633860        B_fluid_int
        62      1       72513120        73992170        B_fluid_int
        70      1       80781443        82123368        B_fluid_int
        155     1       201067953       202583884       B_fluid_int
        230     2       26894103        28819510        B_fluid_int
        249     2       44104111        45189468        B_fluid_int
        266     2       57952946        59251996        B_fluid_int


------------------------------------
[1]: <https://www.nature.com/articles/s41588-022-01017-y>
[2]: <https://github.com/josefin-werme/LAVA?tab=readme-ov-file>
[3]: <https://surfdrive.surf.nl/files/index.php/s/rtmNm8YZERGwl7f>
[4]: <https://surfdrive.surf.nl/s/rtmNm8YZERGwl7f?dir=/cluster_setup/data>
