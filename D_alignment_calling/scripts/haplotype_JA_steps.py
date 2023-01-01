rule plotbams:
    input:  expand(JA_OUT+"bams/{x}.fin.recal.bam.bai", x=samples)
    output: JA_OUT+PBPRE+"plotbams.recalibrated.pdf"
    params: mem="240g",time="8:00:00",partition="ccr",gres="lscratch:200"
    threads: 1
    shell: """
           module load R
           Rscript {PLOT_BAM} {JA_OUT}bams {output} {GENOME}.fasta
           """


rule s2t:
    input:  JA_OUT+"bams/{x}.fin.recal.bam",JA_OUT+"bams/{x}.fin.recal.bam.bai"
    output: temp(JA_OUT+"sam2tsv/{x}.sam2tsv.out")
    params: mem="32g",time="1:00:00",partition="ccr",gres="lscratch:10", sam2tsv=SAM2TSV, genome=GENOME+".fasta"
    threads: 2
    shell: """
           java -jar {params.sam2tsv} -r {params.genome} {input[0]} > {output}
           """

rule t2R:
    input: JA_OUT+"sam2tsv/{x}.sam2tsv.out"
    output: JA_OUT+"sam2tsv/{x}.sam2tsv.pre.RDS"
    threads: 1
    params: mem="20g",time="01:00:00",partition="ccr,norm",gres="lscratch:100"
    shell: """
           module load R
           Rscript {MAKE_PRE_S2T} {input} {output}
           """
           
rule R2p:
    input:  JA_OUT+"sam2tsv/{x}.sam2tsv.pre.RDS"
    output: JA_OUT+"bams/{x}.plotData.CUT.RDS"
    threads: 1
    params: mem="10g",time="1:00:00",partition="ccr,norm",gres="lscratch:100"
    shell: """
           module load R
           Rscript {CUT_HAPL} {input} {output} {SAMPLESHEET} {N_READS_CUT}
           """
                
rule R2plots:
    input: expand(JA_OUT+"bams/{x}.plotData.CUT.RDS", x=samples)
    output: JA_OUT+PBPRE+"_plots.CUT.V2.pdf"
    threads: 1
    params: mem="10g",time="04:00:00",partition="ccr,norm",gres="lscratch:100"
    shell: """
           module load R
           Rscript {MAKE_SPREAD_PLOTS} {JA_OUT}bams/ {SAMPLESHEET} {output}
           """

rule Diffs:
    input:  JA_OUT+"sam2tsv/{x}.sam2tsv.pre.RDS"
    output: JA_OUT+"diffs/{x}.diffs.RDS"
    threads: 1
    params: mem="10g",time="1:00:00",partition="ccr,norm",gres="lscratch:10"
    shell: """
           module load R
           Rscript {DIFFS} {input} {output}
           """