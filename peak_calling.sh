
mkdir -p outputPeaks

bam_dir="/fast/AG_Bunina/Yusuf/Project_Endothelial_and_Stroke/Datasets/HCMEC/atac_output/5.FINAL_OUTPUT"

BAM_FILE="${bam_dir}/D73.final_s.bam"
SAMPLE="D73"

macs2_loc="/fast/AG_Bunina/software_2/macs2_env/bin/python /fast/AG_Bunina/software_2/macs2_env/bin/macs2"

${macs2_loc} callpeak \
  -t "$BAM_FILE" \
  -f BAM -g hs -n "$SAMPLE" \
  --outdir outputPeaks/ \
  --nomodel --shift -100 --extsize 200 -q 0.01
