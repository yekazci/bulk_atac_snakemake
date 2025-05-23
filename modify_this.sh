#!/bin/bash

# User-defined new paths:

# CHANGE THESE PARAMETERS, NOT EXCEED ONE LINE:

main_dir="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac"

# set below true for a dryrun:

dryRun=false

#####################

input_dir="${main_dir}/test_data" 

outdir="${main_dir}/test_snakemake_new"

samples=(D19 D20)                                   # for samples D19_R1.fastq, D19_R2.fastq, D20_R1.fastq, D20_R2.fastq etc.



configFile="${main_dir}/config.yaml"
clusterConfig="${main_dir}/cluster_resource.yaml"    
snakefile="${main_dir}/Snakefile_PE.sh" 
nCores=8

################################### OPTIONAL ##########################

trimmomaticLoc="/fast/AG_Bunina/software_2/Trimmomatic-0.39/trimmomatic-0.39.jar"
trim_adapt="/fast/AG_Bunina/software_2/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
GENOMELOC="/fast/AG_Bunina/annotations/hg38/GRCh38_noalt_as/GRCh38_noalt_as"
tss_bed_file="${main_dir}/human_refs/ensembl_tss_GRCh38.107.bed"
blacklist_bed_file="${main_dir}/human_refs/hg38-blacklist.v2.bed"


############# DO NOT CHANGE BELOW ########################


# --- Override parameters in second_file (params.sh) ---

# File to update
second_file="${main_dir}/params.sh"

# Update dryRun line if it exists; otherwise append it
if grep -q '^dryRun=' "$second_file"; then
  # Replace existing dryRun line
  sed -i "s/^dryRun=.*/dryRun=$dryRun/" "$second_file"
else
  # Append dryRun line if not present
  echo "dryRun=$dryRun" >> "$second_file"
fi

# Update or insert configFile
if grep -q '^configFile=' "$second_file"; then
  sed -i "s|^configFile=.*|configFile=\"$configFile\"|" "$second_file"
else
  echo "configFile=\"$configFile\"" >> "$second_file"
fi

# Update or insert clusterConfig
if grep -q '^clusterConfig=' "$second_file"; then
  sed -i "s|^clusterConfig=.*|clusterConfig=\"$clusterConfig\"|" "$second_file"
else
  echo "clusterConfig=\"$clusterConfig\"" >> "$second_file"
fi

# Update or insert snakefile
if grep -q '^snakefile=' "$second_file"; then
  sed -i "s|^snakefile=.*|snakefile=\"$snakefile\"|" "$second_file"
else
  echo "snakefile=\"$snakefile\"" >> "$second_file"
fi

# Update or insert nCores
if grep -q '^nCores=' "$second_file"; then
  sed -i "s|^nCores=.*|nCores=$nCores|" "$second_file"
else
  echo "nCores=$nCores" >> "$second_file"
fi

###############################################################################

# File to modify
third_file="${main_dir}/runSnakefile.sh"

# Replace the value of params_file in third.sh
if grep -q '^[[:space:]]*params_file=' "$third_file"; then
  sed -i "s|^\([[:space:]]*params_file=\).*|\1\"$second_file\"|" "$third_file"
else
  echo "params_file=\"$second_file\"" >> "$third_file"
fi


# File paths
config_file="${main_dir}/config.yaml"
temp_file=$(mktemp)

in_samples_block=0

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

        # Replace trimmomaticLoc
    if [[ $line =~ ^[[:space:]]*trimmomaticLoc[[:space:]]*:[[:space:]]* ]]; then
        echo "  trimmomaticLoc: \"$trimmomaticLoc\"" >> "$temp_file"
        continue
    fi
    
    # Replace trim_adapt
    if [[ $line =~ ^[[:space:]]*trim_adapt[[:space:]]*:[[:space:]]* ]]; then
        echo "  trim_adapt: \"$trim_adapt\"" >> "$temp_file"
        continue
    fi
    
    # Replace GENOMELOC
    if [[ $line =~ ^[[:space:]]*GENOMELOC[[:space:]]*:[[:space:]]* ]]; then
        echo "  GENOMELOC: \"$GENOMELOC\"" >> "$temp_file"
        continue
    fi
    
    # Replace tss_bed_file
    if [[ $line =~ ^[[:space:]]*tss_bed_file[[:space:]]*:[[:space:]]* ]]; then
        echo "  tss_bed_file: \"$tss_bed_file\"" >> "$temp_file"
        continue
    fi
    
    # Replace blacklist_bed_file
    if [[ $line =~ ^[[:space:]]*blacklist_bed_file[[:space:]]*:[[:space:]]* ]]; then
        echo "  blacklist_bed_file: \"$blacklist_bed_file\"" >> "$temp_file"
        continue
    fi


    # Detect start of samples block
    if [[ $line =~ ^[[:space:]]*samples: ]]; then
        # Write new samples in one-line YAML list
        echo -n "samples: [" >> "$temp_file"
        for ((i=0; i<${#samples[@]}; i++)); do
            if (( i > 0 )); then echo -n ", " >> "$temp_file"; fi
            echo -n "${samples[i]}" >> "$temp_file"
        done
        echo "]" >> "$temp_file"
        echo "" >> "$temp_file"  # Add blank line
        in_samples_block=1
        continue
    fi

    # Skip old samples block lines
    if [[ $in_samples_block -eq 1 ]]; then
        # Stop skipping when we hit a new section or parameter
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

# Replace the original file
mv "$temp_file" "$config_file"

################################################################

sbatch runSnakefile.sh
