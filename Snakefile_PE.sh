#######################################
# General stuff to make things easier #
#######################################
# Make the output nicer and easier to follow
ruleDisplayMessage = "\n\n########################\n# START EXECUTING RULE #\n########################\n"

###########################################
# Onstart, onsuccess and onerror handlers #
###########################################
# Sometimes, it is necessary to specify code that shall be executed when the workflow execution is finished (e.g. cleanup, or notification of the user).
# The onsuccess handler is executed if the workflow finished without error.
onsuccess:
    print("\n\n###############################\n# Workflow finished, no error #\n###############################\n\n")
# Else, the onerror handler is executed.
onerror:
    print("\n\n#####################\n# An error occurred #\n#####################\n\n")
    #shell("mail -s "an error occurred" daria.bunina@mdc-berlin.de < {log}")
# onstart handler will be executed before the workflow starts. Note that dry-runs do not trigger any of the handlers
onstart:
    print("Reading samples and metadata....\n")
    print ("Running workflow for the following samples:\n " + ' \n '.join(map(str, allSamplesUnique)))



#############################
# DIRECTORIES AND VARIABLES #
#############################

#import pandas as pd

# Get all unique sample name:
allSamplesUnique = config["samples"]

# pairedEnd = config["par_general"]["pairedEnd"]

ROOT_dir        = config["par_general"]["outdir"]
INPUT_DIR       = config["par_general"]["input_dir"]
TRIM_dir        = ROOT_dir + "/1.Trimming"
FASTQC_dir      = ROOT_dir + "/2.FastQC"
ALIGN_dir       = ROOT_dir + "/3.ALIGN"
POSTALIGN_dir   = ROOT_dir + "/4.Postalignment"
MAPQ_dir        = POSTALIGN_dir + "/MAPQdir"
chrM_dir        = POSTALIGN_dir + "/chrM"
RSS_dir         = POSTALIGN_dir + "/adjustRSS"
RMdup_dir       = POSTALIGN_dir + "/RMdup"
rmINDEL_dir     = POSTALIGN_dir + "/rmINDEL"
FINAL_dir       = ROOT_dir + "/5.FINAL_OUTPUT"
MACS2_dir       = ROOT_dir + "/6.MACS2_peaks"
POSTALIGN_QC_dir   = ROOT_dir + "/7.Postalignment_QC"
FeatureCounts_dir  = ROOT_dir + "/8.FeatureCounts"
TSS_enrichment_dir = POSTALIGN_QC_dir + "/TSS_results"
TSS_matrix_dir     = TSS_enrichment_dir + "/TSS_matrices"
TSS_plot_dir       = TSS_enrichment_dir + "/TSS_enrichment_plots"
Fragment_dir       = POSTALIGN_QC_dir + "/Frag_length_dist"
FRiP_dir           = POSTALIGN_QC_dir + "/Frip_scores"
Fingerprint_dir    = POSTALIGN_QC_dir + "/Fingerprint_scores"

LOG_BENCHMARK_dir  = ROOT_dir + "/Logs_and_Benchmarks"


#########
# RULES #
#########


rule all:
  input:
    expand("{dir}/{sample}" + "_2.trimmed_fastqc.zip", dir = FASTQC_dir, sample = allSamplesUnique),
    expand("{dir}/{sample}" + ".final_output_bam.bw", dir = FINAL_dir,  sample = allSamplesUnique),
    expand("{dir}/{sample}_peaks.narrowPeak", dir = MACS2_dir, sample = allSamplesUnique),
    MACS2_dir + "/consensus_peaks.bed",
    MACS2_dir + "/consensus_peaks.filtered.bed",
    expand(TSS_matrix_dir + "/{sample}.matrix.gz", sample=allSamplesUnique),
    expand(TSS_plot_dir + "/{sample}_TSS_enrichment.png", sample=allSamplesUnique),
    expand(Fragment_dir + "/{sample}_fragment_length.png", sample=allSamplesUnique),
    expand(FRiP_dir + "/{sample}_frip.txt", sample=allSamplesUnique),
    expand(Fingerprint_dir + "/{sample}_fingerprint_plot.png", sample=allSamplesUnique),
    FRiP_dir + "/frip_scores.tsv",
    FeatureCounts_dir + "/consensus_filtered_peaks.saf",
    FeatureCounts_dir + "/counts.txt"


rule adapter_trim:
  input:
    seqfile=INPUT_DIR + "/{sample}_R1.fastq.gz",
    seqfile2=INPUT_DIR + "/{sample}_R2.fastq.gz"
  output:
    trimmed=TRIM_dir + "/{sample}_1.trimmed.fq",
    trimmed2=TRIM_dir + "/{sample}_2.trimmed.fq"
  params:
    trimmLoc = config["adapters"]["trimmomaticLoc"],
    adapt = config["adapters"]["trim_adapt"],
    illumina = config["par_trimming"]["ILLUMINACLIP_par"],
    trailing = config["par_trimming"]["TRAILING"],
    minlen = config["par_trimming"]["MINLEN"],
    phred = config["par_trimming"]["phred"]
  log:
    LOG_BENCHMARK_dir + "/trimming.{sample}.log"
  message:
        "{ruleDisplayMessage}Trimming of adapters with TRIMMOMATIC in the PE mode for files {input.seqfile:q} and {input.seqfile2:q} using adapters file {params.adapt:q}. This may take a while..."
  threads: 4
  shell:
    """java -jar {params.trimmLoc} PE -threads {threads} -{params.phred} -trimlog {log} {input.seqfile} {input.seqfile2} {output.trimmed} {TRIM_dir}/{wildcards.sample}_R1_unpaired.fq {output.trimmed2} {TRIM_dir}/{wildcards.sample}_R2_unpaired.fq ILLUMINACLIP:{params.adapt:q}:{params.illumina} TRAILING:{params.trailing} MINLEN:{params.minlen}"""

rule fastqc:
  input:
    r1=rules.adapter_trim.output.trimmed,
    r2=rules.adapter_trim.output.trimmed2
  output:
    r1=FASTQC_dir + "/{sample}_1.trimmed_fastqc.zip",
    r2=FASTQC_dir + "/{sample}_2.trimmed_fastqc.zip"
  log:
    LOG_BENCHMARK_dir + "/fastqc.{sample}.log"
  message:
        "{ruleDisplayMessage}Perform FASTQC on the samples {input:q} ..."
  threads: 4
  shell:
    "fastqc -o {FASTQC_dir:q} {input.r1} -t {threads} && fastqc -o {FASTQC_dir:q} {input.r2} -t {threads}"

rule ALIGN:
  input:
    r1=rules.adapter_trim.output.trimmed,
    r2=rules.adapter_trim.output.trimmed2
  output:
    bam = ALIGN_dir + "/{sample}.bt2_s.bam"
  params:
    ref_genome = config["par_align"]["GENOMELOC"],
    window_search = config["par_align"]["WindowSearchPE"],
  log:
    LOG_BENCHMARK_dir + "/alignment.{sample}.log"
  message:
        "{ruleDisplayMessage}Do alignment for files {input:q} with Bowtie2. This may take a while..."
  threads: 8
  shell:
    """
    bowtie2 -p {threads} -X {params.window_search} --very-sensitive -t -x {params.ref_genome} \
    -1 {input.r1} -2 {input.r2} 2> {log} | \
    samtools view -bS - > {output.bam}
    """

#### samtools rules ####

rule samtools_MAPQ_filter:
  input:
    bam = rules.ALIGN.output.bam
  output:
    MAPQsort = MAPQ_dir + "/{sample}.bt2_sort.bam",
    index = MAPQ_dir + "/{sample}.bt2_sort.bam.bai",
    MAPQf = MAPQ_dir + "/{sample}.q10.bt2_f.bam",
    index2 = MAPQ_dir + "/{sample}.q10.bt2_f.bam.bai",
    stats = MAPQ_dir + "/{sample}.stats",
    csv   = MAPQ_dir + "/{sample}.csv.gz"
  params:
    MAPQ = config["par_postalign"]["MAPQ"]
  log:
    LOG_BENCHMARK_dir + "/MAPQf.{sample}.log"
  threads: 4
  shell:
    """
    samtools sort -o {output.MAPQsort:q} --threads {threads} {input.bam:q} &&
    samtools index {output.MAPQsort:q}  &&
    samtools view -b -q {params.MAPQ} -F 4 {output.MAPQsort:q} > {output.MAPQf:q}  &&
    samtools index {output.MAPQf:q}  &&
    samtools flagstat {output.MAPQf:q} > {output.stats:q}  &&
    samtools view {output.MAPQf:q} | cut -f3,5,9 | gzip  -f > {output.csv:q}
    """

rule samtools_remove_chrM:
  input:
    bam = rules.samtools_MAPQ_filter.output.MAPQf
  output:
    final_bam = chrM_dir + "/{sample}.q10.rmcM.bt2_s.bam",
    index2 = chrM_dir + "/{sample}.q10.rmcM.bt2_s.bam.bai",
    stats = chrM_dir + "/{sample}.q10.rmcM.bt2_s.bam.stats"
  threads: 1
  shell:
    """
    samtools idxstats {input.bam:q} | cut -f 1 | grep chr | grep -Pv "chrM|chrUn|random|hap" | xargs samtools view -b {input.bam:q} > {output.final_bam:q} &&
    samtools index {output.final_bam:q} &&
    samtools flagstat {output.final_bam:q} > {output.stats:q}
    """

rule remove_duplicates:
  input:
    adRSS_bam = rules.samtools_remove_chrM.output.final_bam
  output:
    rmDUP_bam = RMdup_dir + "/{sample}.q10.rmcM.rmDUP.bt2_s.bam",
    flagstat_txt = RMdup_dir + "/{sample}.ST.rmDUP.txt",
    index = RMdup_dir + "/{sample}.q10.rmcM.rmDUP.bt2_s.bam.bai"
  log:
    rmDup_log = RMdup_dir + "/{sample}.rmDup.log"
  threads: 1
  shell:
    """
    samtools rmdup -s {input.adRSS_bam:q} {output.rmDUP_bam} && samtools index {output.rmDUP_bam} &&
    samtools flagstat {output.rmDUP_bam} > {output.flagstat_txt}
    """

# the adjust_RSS rule does perform shifting of reads to account for the Tn5 transposase's 9 bp binding offset.

# This +4/-5 adjustment is a standard approximation to correct for Tn5's 9 bp binding footprint.

rule adjust_RSS:
  input:
    rmcM_bam = rules.remove_duplicates.output.rmDUP_bam
  output:
    bam = FINAL_dir + "/{sample}.final.bam",
    stats = RSS_dir + "/{sample}.ST.adRSS.txt",
    csv = RSS_dir + "/{sample}.adRSS.csv"
  threads: 1
  shell:
    """
    time \
    cat <(samtools view -H {input.rmcM_bam:q}) <(samtools view -F 16 {input.rmcM_bam:q} | awk 'BEGIN {{OFS = "\\t"}} ; {{$4=$4+4; print $0}}') <(samtools view -f 16 {input.rmcM_bam:q} | awk 'BEGIN {{OFS = "\\t"}} ; {{$4=$4-5; print $0}}') | samtools view -S -b -o {output.bam:q} - &&
        samtools flagstat {output.bam:q} > {output.stats:q} &&
        samtools view {output.bam:q} | cut -f3,5,9 | gzip  -f > {output.csv:q}
    """

# ChatGPT note: 

# There is a commented out picard_sort rule (Picard SortSam) — probably initially considered but ultimately replaced by samtools sort.

# rule picard_sort:
#   input:
#     rules.adjust_RSS.output.bam
#   output:
#     bam_final = FINAL_dir + "/{sample}.final_s.bam",
#     index = FINAL_dir + "/{sample}.final_s.bam.bai"
#   params:
#     picardLoc = config["par_postalign"]["PicardLoc"]
#   message:
#     "{ruleDisplayMessage}Sort bam file {input:q} with PicardTools..."
#   log:
#     LOG_BENCHMARK_dir + "/{sample}.sortPicard.log"
#   threads: 1
#   shell:
#     """
#     java -jar {params.picardLoc} SortSam \
#                     I={input:q} \
#                     O={output.bam_final:q}  \
#                     SORT_ORDER=coordinate \
#                     VALIDATION_STRINGENCY=SILENT \
#                     2> {log:q} &&
#     samtools index {output.bam_final}
#     """

rule bigwig:
  input:
    rules.adjust_RSS.output.bam
  output:
    bam_final = FINAL_dir + "/{sample}.final_s.bam",
    index = FINAL_dir + "/{sample}.final_s.bam.bai",
    bw = FINAL_dir + "/{sample}.final_output_bam.bw"
  message:
    "{ruleDisplayMessage}Convert bam file {input:q} to bigwig..."
  threads: 4
  shell:
    """
    samtools sort -o {output.bam_final:q} --threads {threads} {input:q} &&
    samtools index {output.bam_final:q} &&
    bamCoverage -b {output.bam_final:q} -o {output.bw} --binSize 50 --effectiveGenomeSize 2864785220 --normalizeUsing RPGC --numberOfProcessors {threads} --minMappingQuality 10
    """

rule call_peaks_macs2:
  input:
    bam = FINAL_dir + "/{sample}.final_s.bam"
  params:
    sampleName = "{sample}",
    macs2 = config["par_postalign_qc"]["macs2_loc"]
  output:
    peak = MACS2_dir + "/{sample}_peaks.narrowPeak"
  log:
    out = LOG_BENCHMARK_dir + "/macs2.{sample}.out",
    err = LOG_BENCHMARK_dir + "/macs2.{sample}.err"
  message:
    "{ruleDisplayMessage}Calling peaks on {input.bam:q} with MACS2..."
  threads: 1
  shell:
    """
    {params.macs2} callpeak \
      -t {input.bam} \
      -f BAM -g hs -n {params.sampleName} \
      --outdir {ROOT_dir}/6.MACS2_peaks \
      --nomodel --shift -100 --extsize 200 -q 0.01 \
      > {log.out} 2> {log.err}
    """

rule consensus_peaks:
  input:
    expand(MACS2_dir + "/{sample}_peaks.narrowPeak", sample=allSamplesUnique)
  output:
    peaks = MACS2_dir + "/consensus_peaks.bed",
    filtered_peaks = MACS2_dir + "/consensus_peaks.filtered.bed"
  params:
    blacklist = config["par_postalign_qc"]["blacklist_bed_file"]
  log:
    LOG_BENCHMARK_dir + "/consensus_peaks.log"
  message:
    "{ruleDisplayMessage}Merging all MACS2 narrowPeak files to create consensus peaks..."
  threads: 1
  shell:
    """
    cat {input} | sort -k1,1 -k2,2n | bedtools merge -i - > {output.peaks}

    bedtools intersect -v -a {output.peaks} -b {params.blacklist} > {output.filtered_peaks}
    """

rule tss_enrichment:
  input:
    bw = FINAL_dir + "/{sample}.final_output_bam.bw"
  params:
    tss = config["par_postalign_qc"]["tss_bed_file"]
  output:
    matrix = TSS_matrix_dir + "/{sample}.matrix.gz",
    plot = TSS_plot_dir + "/{sample}_TSS_enrichment.png"
  threads: 4
  log:
    LOG_BENCHMARK_dir + "/TSS_enrichment.{sample}.log"
  message:
    "{ruleDisplayMessage}Compute TSS enrichment matrix and plot for {input.bw:q} ..."
  shell:
    """
    computeMatrix reference-point \
      --referencePoint TSS \
      -b 2000 -a 2000 \
      -R {params.tss} \
      -S {input.bw} \
      --skipZeros \
      --sortRegions descend \
      --sortUsing mean \
      -o {output.matrix} \
      -p {threads} > {log} 2>&1

    plotProfile -m {output.matrix} -out {output.plot}
    """

rule fingerprint_plot:
  input:
    bam = FINAL_dir + "/{sample}.final_s.bam",
  output:
    plot = Fingerprint_dir + "/{sample}_fingerprint_plot.png"
  threads: 1
  log:
    LOG_BENCHMARK_dir + "/{sample}_Fingerprint_plot.log"
  message:
    "{ruleDisplayMessage}Generate fingerprint plot from final BAMs ..."
  shell:
    """
    plotFingerprint -b {input.bam} \
    -plot {output.plot} > {log} 2>&1
    """


rule fragment_length_distribution:
  input:
    bam = FINAL_dir + "/{sample}.final_s.bam"
  params:
    sampleName = "{sample}"
  output:
    png = Fragment_dir + "/{sample}_fragment_length.png"
  threads: 4
  log:
    LOG_BENCHMARK_dir + "/Fragment_length.{sample}.log"
  message:
    "{ruleDisplayMessage}Compute fragment length distribution for {input.bam:q} ..."
  shell:
    """
    bamPEFragmentSize -b {input.bam} \
      --verbose \
      --maxFragmentLength 1000 \
      --histogram {output.png} \
      --plotTitle {params.sampleName} \
      --samplesLabel "" \
      --binSize 500 \
      --distanceBetweenBins 500000 \
      -p {threads} > {log} 2>&1
    """

rule frip_score:
  input:
    bam = FINAL_dir + "/{sample}.final_s.bam",
    peaks = MACS2_dir + "/consensus_peaks.filtered.bed"
  output:
    score = FRiP_dir +  "/{sample}_frip.txt"
  log:
    LOG_BENCHMARK_dir + "/FRiP_score.{sample}.log"
  message:
    "{ruleDisplayMessage}Calculate FRiP score for {input.bam:q} ..."
  shell:
    """
    FRIP=$(echo "scale=4; \
    $(bedtools intersect -a "{input.bam}" -b "{input.peaks}" -bed | wc -l) / \
    $(samtools view -c "{input.bam}")" | bc)

    echo "$FRIP" > "{output.score}"
    """

rule merge_frip_scores:
  input:
    expand(FRiP_dir + "/{sample}_frip.txt", sample=config["samples"])
  output:
    FRiP_dir + "/frip_scores.tsv"
  message:
    "{ruleDisplayMessage}Merge all individual FRiP scores into a single file ..."
  shell:
    """
    for f in {FRiP_dir}/D*_frip.txt; do 
        sample=$(basename "$f" _frip.txt); 
        score=$(grep -o '[0-9.]*$' "$f"); 
        echo -e "${{sample}}\t${{score}}"; 
    done | sort > {output}
    """

rule featureCounts:
  input:
    bams = expand(FINAL_dir + "/{sample}.final_s.bam", sample=config["samples"]),
    peaks = MACS2_dir + "/consensus_peaks.filtered.bed"
  output:
    saf = FeatureCounts_dir + "/consensus_filtered_peaks.saf",
    counts = FeatureCounts_dir + "/counts.txt"
  log:
    LOG_BENCHMARK_dir + "/featureCounts.log"
  threads: 16
  message:
    "{ruleDisplayMessage}Generate count matrix..."
  shell:
    """
    awk 'BEGIN{{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"}} \
    {{print "peak"NR, $1, $2, $3, "."}}' \
    {input.peaks} > {output.saf}

     featureCounts -F SAF \
      -a {output.saf} \
      -o {output.counts} \
      -T {threads} \
      --read2pos 5 \
      -p \
      {input.bams} > {log} 2>&1
    
    """

      

