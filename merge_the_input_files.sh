
# One can use the following code to merge the fastq files from different lanes.
# It also simplifies their naming format for easy addressing in the snakemake file.

# Example raw files:

# 5055_D_01_S22_L006_R1_001.fastq.gz

# A5055_D_01_S22_L006_R2_001.fastq.gz

# A5055_D_01_S22_L007_R1_001.fastq.gz

# A5055_D_01_S22_L007_R2_001.fastq.gz

# Now, for example, D_01.....R1 files are merged into D01_R1.fastq.gz

# Below, this merging is done for all R1 and R2 reads of all samples.

# You can adjust it for your naming format.

# Create output directory if it doesn't exist
mkdir -p merged_fastq

# example 

# Loop from 01 to 96
for i in $(seq -w 1 96); do
    sample_id="D_$i"
    merged_prefix="D${i}"

    echo "Merging files for sample: $sample_id"

    cat $(ls A5055_${sample_id}*_R1_001.fastq.gz | sort) > merged_fastq/${merged_prefix}_R1.fastq.gz
    cat $(ls A5055_${sample_id}*_R2_001.fastq.gz | sort) > merged_fastq/${merged_prefix}_R2.fastq.gz
done
