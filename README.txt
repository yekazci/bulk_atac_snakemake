
User should only modify the "modify_this.sh" file.

There are five variables to set in this file:

# set the following variable to the folder path that contains pipeline scripts:

main_dir="/fast/AG_Bunina/Yusuf/ATAC_BULK_snakemake/bunina_lab_atac"

# set below true for a dryrun, highly recommended:

dryRun=false

# set below the input folder path.

# If you copy links of your fastq files to the test_data folder, you can skip this.

input_dir="${main_dir}/test_data" 

# OPTIONAL: set below the output folder path. 

outdir="${main_dir}/test_snakemake_new"

# IMPORTANT:

# For samples that are named as below:

D01_R1.fastq
D01_R2.fastq
D02_R1.fastq
D02_R2.fastq

# Write the identifier prefixes for each sample like below:

samples=(D01 D02)