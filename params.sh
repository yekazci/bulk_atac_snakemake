# Config file to use
### ADJUST THE PATH TO YOUR CONFIG FILE LOCATION!###
# STRING. Absolute path to the main pipeline configuration YAML file.

# CHANGE HERE to the directory that contains the pipeline scripts.

# configFile="/fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/HCMEC/bunina_lab_atac/config.yaml"

configFile="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac/config.yaml"

# Maximum number of CPUs per rule executed. See the separate cluster specification
# INTEGER. Maximum number of cores (threads) to be used by Snakemake per rule.
# nCores=8 # set this to 32 for large sample size.

nCores=8

# Use a dry run for testing purposes?
# BOOLEAN. If set to "true", Snakemake will simulate the workflow without executing any jobs.
dryRun=true

### Sections below: DO NOT CHANGE unless you know what you do ###

# Snakefile to use
# STRING. Path to the Snakefile that contains the workflow logic for paired-end ATAC-seq analysis.

# CHANGE HERE to the directory that contains the pipeline scripts.

# snakefile="/fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/HCMEC/bunina_lab_atac/Snakefile_PE.sh"
snakefile="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac/Snakefile_PE.sh"

# Cluster options:
# YEK : The following purports to provide path to a json file, but the file is in yaml format. So the extension is wrongly written ".json". I corrected it below.
# STRING. Path to the JSON file specifying cluster resource usage for each Snakemake rule (e.g., memory, CPUs, runtime).

# CHANGE HERE to the directory that contains the pipeline scripts.

# clusterConfig="/fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/HCMEC/bunina_lab_atac/cluster_resource.yaml"        
clusterConfig="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac/cluster_resource.yaml"