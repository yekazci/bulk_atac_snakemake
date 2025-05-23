#!/bin/bash

# User-defined new paths:

# CHANGE THESE PARAMETERS, NOT EXCEED ONE LINE:

main_dir="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac_test"  # Main working directory for the ATAC-seq pipeline

# set below true for a dryrun:

dryRun=true                                                           # Set to true to perform a Snakemake dry-run without executing jobs

# Input folder:

input_dir="${main_dir}/test_data"  # Directory containing the input FASTQ files

# Input names:

samples=(D19 D20)  # Sample identifiers used in the pipeline (e.g., D19_R1.fastq)

################################### OPTIONAL: CHANGE CAREFULLY ##########################

# output_folder

outdir="${main_dir}/test_snakemake_new"  # Output directory for results

# Optional tool and reference paths
trimmomaticLoc="/fast/AG_Bunina/software_2/Trimmomatic-0.39/trimmomatic-0.39.jar"  # Trimmomatic JAR file path
trim_adapt="/fast/AG_Bunina/software_2/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"  # Adapter sequences for trimming
GENOMELOC="/fast/AG_Bunina/annotations/hg38/GRCh38_noalt_as/GRCh38_noalt_as"  # Reference genome fasta basename
tss_bed_file="${main_dir}/human_refs/ensembl_tss_GRCh38.107.bed"  # BED file with transcription start sites
blacklist_bed_file="${main_dir}/human_refs/hg38-blacklist.v2.bed"  # BED file with blacklisted regions

nCores=8  # Number of CPU cores to allocate for Snakemake execution

############# DO NOT CHANGE BELOW ########################

# Configuration files and Snakemake parameters
configFile="${main_dir}/config.yaml"  # Snakemake configuration YAML file
clusterConfig="${main_dir}/cluster_resource.yaml"  # Cluster resources YAML for job submission
snakefile="${main_dir}/Snakefile_PE.sh"  # The main Snakemake workflow file

# --- Override parameters in second_file (params.sh) ---

second_file="${main_dir}/params.sh"  # Parameter shell script file to be updated

# Update or append dryRun variable
if grep -q '^dryRun=' "$second_file"; then
  sed -i "s/^dryRun=.*/dryRun=$dryRun/" "$second_file"  # Replace existing value
else
  echo "dryRun=$dryRun" >> "$second_file"  # Append if not present
fi

# Update or insert configFile path
if grep -q '^configFile=' "$second_file"; then
  sed -i "s|^configFile=.*|configFile=\"$configFile\"|" "$second_file"
else
  echo "configFile=\"$configFile\"" >> "$second_file"
fi

# Update or insert clusterConfig path
if grep -q '^clusterConfig=' "$second_file"; then
  sed -i "s|^clusterConfig=.*|clusterConfig=\"$clusterConfig\"|" "$second_file"
else
  echo "clusterConfig=\"$clusterConfig\"" >> "$second_file"
fi

# Update or insert snakefile path
if grep -q '^snakefile=' "$second_file"; then
  sed -i "s|^snakefile=.*|snakefile=\"$snakefile\"|" "$second_file"
else
  echo "snakefile=\"$snakefile\"" >> "$second_file"
fi

# Update or insert nCores value
if grep -q '^nCores=' "$second_file"; then
  sed -i "s|^nCores=.*|nCores=$nCores|" "$second_file"
else
  echo "nCores=$nCores" >> "$second_file"
fi

###############################################################################

# File to modify
third_file="${main_dir}/runSnakefile.sh"  # Snakemake run script to be updated

# Replace the value of params_file in third.sh
if grep -q '^[[:space:]]*params_file=' "$third_file"; then
  sed -i "s|^\([[:space:]]*params_file=\).*|\1\"$second_file\"|" "$third_file"
else
  echo "params_file=\"$second_file\"" >> "$third_file"
fi


############## use same logic to assign wrapper_file to the fourth_file path:

fourth_file="${main_dir}/runSnakemakeWrapper.sh"  # Wrapper script whose path we want to assign

# Replace the value of wrapper_file in runSnakefile.sh
if grep -q '^[[:space:]]*wrapper_file=' "$third_file"; then
  sed -i "s|^\([[:space:]]*wrapper_file=\).*|\1\"$fourth_file\"|" "$third_file"
else
  echo "wrapper_file=\"$fourth_file\"" >> "$third_file"
fi


# Modify configuration YAML to reflect new input/output paths and references
config_file="${main_dir}/config.yaml"
temp_file=$(mktemp)  # Create a temporary file to hold updated config

in_samples_block=0  # Flag to track if we are inside the 'samples' YAML block

while IFS= read -r line; do
    # Replace input_dir
    if [[ $line =~ ^[[:space:]]*input_dir: ]]; then
        echo "  input_dir: \"$input_dir\"" >> "$temp_file"
        continue
    fi

    # Replace outdir
    if [[ $line =~ ^[[:space:]]*outdir: ]]; then
        echo "  outdir: \"$outdir\"" >> "$temp_file"
        continue
    fi

    # Replace tool/reference paths
    if [[ $line =~ ^[[:space:]]*trimmomaticLoc[[:space:]]*:[[:space:]]* ]]; then
        echo "  trimmomaticLoc: \"$trimmomaticLoc\"" >> "$temp_file"
        continue
    fi
    if [[ $line =~ ^[[:space:]]*trim_adapt[[:space:]]*:[[:space:]]* ]]; then
        echo "  trim_adapt: \"$trim_adapt\"" >> "$temp_file"
        continue
    fi
    if [[ $line =~ ^[[:space:]]*GENOMELOC[[:space:]]*:[[:space:]]* ]]; then
        echo "  GENOMELOC: \"$GENOMELOC\"" >> "$temp_file"
        continue
    fi
    if [[ $line =~ ^[[:space:]]*tss_bed_file[[:space:]]*:[[:space:]]* ]]; then
        echo "  tss_bed_file: \"$tss_bed_file\"" >> "$temp_file"
        continue
    fi
    if [[ $line =~ ^[[:space:]]*blacklist_bed_file[[:space:]]*:[[:space:]]* ]]; then
        echo "  blacklist_bed_file: \"$blacklist_bed_file\"" >> "$temp_file"
        continue
    fi

    # Detect start of samples list
    if [[ $line =~ ^[[:space:]]*samples: ]]; then
        echo -n "samples: [" >> "$temp_file"
        for ((i=0; i<${#samples[@]}; i++)); do
            if (( i > 0 )); then echo -n ", " >> "$temp_file"; fi
            echo -n "${samples[i]}" >> "$temp_file"
        done
        echo "]" >> "$temp_file"
        echo "" >> "$temp_file"
        in_samples_block=1
        continue
    fi

    # Skip old samples block until a new key is detected
    if [[ $in_samples_block -eq 1 ]]; then
        if [[ $line =~ ^[[:space:]]*[-#] ]]; then
            continue
        elif [[ $line =~ ^[[:space:]]*$ ]]; then
            continue
        elif [[ $line =~ ^[[:alnum:]_]+: ]]; then
            in_samples_block=0
            echo "$line" >> "$temp_file"
        else
            continue
        fi
    else
        echo "$line" >> "$temp_file"
    fi
done < "$config_file"

# Replace original YAML with the modified version
mv "$temp_file" "$config_file"

################################################################

# Submit the Snakemake job to the scheduler
sbatch runSnakefile.sh
