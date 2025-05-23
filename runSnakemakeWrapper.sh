#!/bin/bash


for dir in \
adapter_trim \
fastqc \
ALIGN \
samtools_MAPQ_filter \
samtools_remove_chrM \
remove_duplicates \
adjust_RSS \
bigwig \
call_peaks_macs2 \
consensus_peaks \
tss_enrichment \
fingerprint_plot \
fragment_length_distribution \
frip_score \
merge_frip_scores \
featureCounts
do
  mkdir -p "snakemake_logs/$dir"
done

# the snakemake pipeline is only working with this snakemake version !!!!!!!!!!!!!!!!!!!!!!!! 

# snakemake_loc="/gnu/store/c2i04qdzah673yk300wv87xalgbms4iw-snakemake-7.7.0/bin"

# the snakemake pipeline is only working with this snakemake version !!!!!!!!!!!!!!!!!!!!!!!! 

snakemake_dir=/fast/AG_Bunina/software_2/snakemake_7_7_0_env

snakemake_loc="${snakemake_dir}/bin/python ${snakemake_dir}/bin/snakemake"

mkdir -p /fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/HCMEC/bunina_lab_atac/test_snakemake_new

# Check if dryRun is set to true
if [ "$dryRun" = true ] ; then
        # If dryRun is true, execute Snakemake with a dry run option to simulate the workflow without running the actual tasks
        ${snakemake_loc} --cluster sbatch --jobs $nCores -s $snakefile --configfile $configFile --dryrun --cores 1
else
        # If dryRun is not true, execute Snakemake with the actual parameters

            # Specify the cluster configuration for Slurm jobs
        ${snakemake_loc} \
    --cluster 'sbatch --cpus-per-task={cluster.nCPUs} --mem={cluster.memory} --time={cluster.maxTime} --job-name={cluster.name} --output={cluster.output} -e {cluster.error}' \
    --jobs 96 \
    -s $snakefile --restart-times 3 \
    --cluster-config $clusterConfig \
    --configfile $configFile \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    
fi



