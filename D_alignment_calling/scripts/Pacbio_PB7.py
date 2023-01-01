import os
### sbatch ../../../../scripts/scripts_batch/batch_2j_ccr.sh PacBio_PB3D.py
### snakemake -s genotype.py --unlock
### sbatch batch.sh genotype.py  ##snakemake on cluster node
### sbatch genotype.py  ## snakemake running on local node

PBPRE = "PB7_align"
SAMPLESHEET = "../B_barcode_report/output/PB7_NH2700_JA_meta.xlsx"
DATA = "../PB7_NH2700_JA_meta_binned_fastqs/"

### JA Custom Scripts and helpers

MAKE_PRE_S2T = "D_alignment_calling/scripts/JA_make_pre.sam2tsv.R"  
CUT_HAPL = "D_alignment_calling/scripts/JA_create.haplotypes.CUT.R"  
MAKE_SPREAD_PLOTS = "D_alignment_calling/scripts/JA_make.spread.plots.V3.R"  
PLOT_BAM = "D_alignment_calling/scripts/JA_plotbam.recalBams.R"  


### N_READS_CUT removes haplotypes less than N and also only top N haplotypes are plotted
N_READS_CUT = 7

GENOME = "D_alignment_calling/scripts/Pacbio_STE50toFUS1/STE50toFUS1_S4921" 

## PACKAGES ON BUIWULF - HARDCODED INPUTS


GATK="java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar /usr/local/apps/GATK/3.5-0/GenomeAnalysisTK.jar"
PICARD="module load picard/2.1.1; java -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups"
SAM2TSV = "/home/ahujajs/modulefiles/jvarkit/dist/sam2tsv.jar"



TMPR = "temp"
FOLDER = PBPRE+"/"
JA_OUT=FOLDER
samples, = glob_wildcards(DATA+"{sample}.fastq")

print("Samples are "+str(samples))
print("samplesheet is "+str(SAMPLESHEET))


rule all_wgslow:
    input: expand(JA_OUT+"bams/{s}.fin.recal.bam.bai",s=samples),
           GENOME+".dict",
#           expand(JA_OUT+"vcfs/{s}.g.vcf",s=samples),
           expand(JA_OUT+"sam2tsv/{x}.sam2tsv.pre.RDS", x=samples),
           expand(JA_OUT+"bams/{x}.plotData.CUT.RDS", x=samples),
           JA_OUT+PBPRE+"_plots.CUT.V2.pdf",
#           expand(JA_OUT+"plotbamraster/{x}.png", x=samples),
           JA_OUT+PBPRE+"plotbams.recalibrated.pdf"
    output: 
    params: rname="final"


include: "D_alignment_calling/scripts/Pacbio_alignment_dedup_ETC_fastq.py"
include: "D_alignment_calling/scripts/haplotype_JA_steps.py"

