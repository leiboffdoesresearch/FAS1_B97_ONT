#load info from config file
#config file contains sample name info, genome file info, and parameters
configfile: "config.yaml"

#use minimal debian containerized environment with conda
#useful for OS standardization
container: "docker://continuumio/miniconda3:4.5.11"

#this rule looks for all the final files
#drives the back propagation of all intermediate rules
#elimintates the need for specifying inputs or outputs on command line
rule all:
    input:
        expand('sniffles/{sample}_{genome}_DP{sniffles_depth}.vcf', 
            sample=config["samples"], 
            genome=config["genomes"], 
            sniffles_depth=config["sniffles_depth"])


# #nanoplot ONT run
# rule nanoplot_ONT:
#     input:
#         R = '{sample}/sequencing_summary/*_sequencing_summary.txt'
#     output:
#         R1 = 'trimmed_reads/{sample}_R1.fastq.gz',
#         R2 = 'trimmed_reads/{sample}_R2.fastq.gz',
#         html_report = 'fastp_report/{sample}.html',
#         json_report = 'fastp_report/{sample}.json'
#     conda:
#         "envs/nanoplot.yaml"
#     threads: 4 #how many?
#     shell: 
#         "NanoPlot -t {threads} --verbose "
#         "-o ./nanoplot/{sample} "
#         "--plots kde hex dot pauvre "
#         "--N50 --title "WAB2 RAD004 ONT 02/26" " 
#         "--summary ./ONT_stats/NanoporT_20190226_FAK63765_MN27201_sequencing_run_WAB2_MUT_Pool_97948_sequencing_summary.txt"

#concatenate all fastq files from ONT output
#future pipes should also use guppy to call reads from raw wiggles
#rule cat_ONT:
#    input:


# #nanoplot reads
# rule nanoplot_fastq:
#     input:
#         R = 'cat_reads/{sample}.fastq.gz'
#     output:
#         R1 = 'trimmed_reads/{sample}_R1.fastq.gz',
#         R2 = 'trimmed_reads/{sample}_R2.fastq.gz',
#         html_report = 'fastp_report/{sample}.html',
#         json_report = 'fastp_report/{sample}.json'
#     conda:
#         "envs/fastp.yaml"
#     threads: 4 #how many?
#     shell: 
#         "fastp -w {threads} "
#         "-i {input.R1} -I {input.R2} "
#         "-o {output.R1} -O {output.R2} "
#         "-h {output.html_report} -j {output.json_report}"

#minimap2 in ONT mode
rule minimap2:
    input:
        R = 'cat_reads/{sample}.fastq'
        genome = 'genomes/{genome}.fa'
    output:
        bam = 'aligned_reads/{sample}_{genome}.bam'
    log: 'aligned_reads/{sample}_{genome}.log'
    conda:
        "envs/minimap2.yaml"
    threads: 20 # how many?
    shell:
        "minimap2 -ax map-ont {input.genome} "
        "-I 10G -t {threads} -2 -K 10G "
        "{input.R} "
        "2> {log} "
        "| samtools view -bS - "
        "> {output.bam}"

# #nanoplot ONT aligned
# rule nanoplot_bam:
#     input:
#         R = 'concat_reads/{sample}.fastq.gz'
#     output:
#         R1 = 'trimmed_reads/{sample}_R1.fastq.gz',
#         R2 = 'trimmed_reads/{sample}_R2.fastq.gz',
#         html_report = 'fastp_report/{sample}.html',
#         json_report = 'fastp_report/{sample}.json'
#     conda:
#         "envs/nanoplot.yaml"
#     threads: 4 #how many?
#     shell: 
#         "fastp -w {threads} "
#         "-i {input.R1} -I {input.R2} "
#         "-o {output.R1} -O {output.R2} "
#         "-h {output.html_report} -j {output.json_report}"


#sort BAMs
rule BAM_sort:
    input:
        bam = 'aligned_reads/{sample}_{genome}.bam'
    output:
        sorted_bam = 'bam_sorted/{sample}_{genome}.sorted.bam'
    conda:
        "envs/samtools.yaml"
    threads: 20
    shell:
        "samtools sort -@ {threads} {input.bam} -o {output.sorted_bam}"

#index BAMs
rule BAM_index:
    input:
        sorted_bam = 'bam_sorted/{sample}_{genome}.sorted.bam'
    output:
        bam_index = 'bam_sorted/{sample}_{genome}.sorted.bam.bai'
    conda:
        "envs/samtools.yaml"
    threads: 20
    shell:
        "samtools index -@ {threads} {input.sorted_bam}"

# create MD tag of IN/DELs using SAMtools
# suppress STDout -- will announce every read
rule samtools_calmd:
    input:
        sorted_bam = 'bam_sorted/{sample}_{genome}.sorted.bam'
        bam_index = 'bam_sorted/{sample}_{genome}.sorted.bam.bai'
        genome = 'genomes/{genome}.fa'
    output:
        MD_bam = 'bam_MD/{sample}_{genome}.MD.bam'
    conda:
        "envs/samtools.yaml"
    threads: 20
    shell:
        "samtools calmd -b {input.sorted_bam} {genome} > {output.MD_bam}"

# run sniffles to detect SVs
# set depth req, s to match sequencing depth
rule sniffles_SV:
    input:
        MD_bam = 'bam_MD/{sample}_{genome}.MD.bam'
        depth = '{depth}'
    output:
        vcf = 'sniffles/{sample}_{genome}_DP{sniffles_depth}.vcf'
    conda:
        "envs/sniffles.yaml"
    threads: 40
    shell:
        "sniffles -s {input.depth} --genotype --report_seq "
        "-n {threads} -m {input.MD_bam} -v {output.vcf}"

##### CUT AND PASTED THOUGHTS ######
# Run SURVIVOR merge files 1000 1 1 0 0 0 output_merged.vcf

# Then you have a multi sample vcf file for all samples. Feel free to alter the parameters (e.g. 1000).
# If you want to take additional steps to be more precise when comparing things (beta version):

# run Sniffles again: sniffles -m reads.bam -v output_genotyped.vcf --Ivcf output_merged.vcf

# This genotypes all the variations that you found previously.
# merge the so generated files again with:

# ls *genotyped.vcf > files_genotyped
# SURVIVOR merge 1000 0 1 0 0 0 output_genotyped_merged.vcf
