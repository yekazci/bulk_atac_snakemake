#!/bin/bash

#SBATCH --job-name=bulk_atac        # Name of the SLURM job
#SBATCH --mail-user=yusufenes.kazci@mdc-berlin.de
#SBATCH --output=slurm_output_%j.txt
#SBATCH --error=slurm_error_%j.txt


#load guix profile

export GUIX_PROFILE="/fast/AG_Bunina/common_guix_profiles/bulk_atac"          # Set the path to the user's Guix profile
source "$GUIX_PROFILE/etc/profile"                                            # Source the environment settings from the Guix profile
export GUIX_LOCPATH="$GUIX_PROFILE/lib/locale"                                # Set the locale path for language and region settings

########################
# PATHS AND PARAMETERS #
########################

# CHANGE HERE to the directory that contains the pipeline scripts.

# All parameters and paths are defined here (UPDATE TO YOUR PATH!!):

params_file="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac_test/params.sh"

. ${params_file} 

# Source external parameter file which defines variables like $nCores, $snakefile, and $configFile

########################
# RUN AUTOMATED SCRIPT #
########################

# snakemake --cluster "sbatch {threads}" --jobs $nCores -s $snakefile $configFile --latency-wait 120
# Example of how Snakemake can be manually run with cluster support; commented out here

# CHANGE HERE to the directory that contains the pipeline scripts.

# Run a wrapper script that internally executes Snakemake with the loaded parameters

wrapper_file="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac_test/runSnakemakeWrapper.sh"

. ${wrapper_file} 

# . "/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac/runSnakemakeWrapper.sh"  

