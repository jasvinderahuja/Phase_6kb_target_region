### jasvinder.ahuja@nih.gov
### jasvinderahuja@gmail.com
### 

rule flagstat:
    input: JA_OUT+"bams/{x}.bam"
    output: JA_OUT+"bams/{x}.bam.flags"
    params:  mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
    threads: 1
    shell: """
            module load samtools
            samtools flagstat {input}>{output}
            """
            
            
rule minimap2:
    input: DATA+"{sample}.fastq"
    output: temp(JA_OUT+"bams/{sample}.minimap2.bam")
    group: "align"
#    threads: 4
    params:  genome=GENOME+".fasta"
    shell: """ 
           module load minimap2
           module load samtools
           minimap2 -ax map-pb --secondary=no {GENOME}.fasta {input} | samtools view -bS - > {output}
           """ 


rule samtools_sort:
    input: JA_OUT+"bams/{x}.minimap2.bam"
    output: temp(JA_OUT+"bams/{x}.sorted.bam")
    group: "align"
#    threads: 1
#    params: rname="sort", mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
    shell: """
    	   module load samtools
           samtools sort -T {wildcards.x} -O bam {input} > {output}
           """


rule faidx:
    input: "{GENOME}.fasta"
    output: "{GENOME}.fasta.fai"
    threads: 1
    params:  mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
    shell: """
           module load samtools
           samtools faidx {input}
           """

rule indexbam2:
    input: JA_OUT+"bams/{x}.fin.recal.bam"
    output: JA_OUT+"bams/{x}.fin.recal.bam.bai"
    group: "align"
#    threads: 1
#    params: mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
    shell: """
           module load samtools
           samtools index {input}
           if [ ! -f {GENOME}.fasta.fai ]; then
               samtools faidx {GENOME}.fasta
              fi
           """
      
rule indexbam:
    input: JA_OUT+"bams/{x}.fin.bam"
    output: temp(JA_OUT+"bams/{x}.fin.bam.bai")
#    threads: 1
	group: "align"
#    params: mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
    shell: """
           module load samtools
           samtools index {input}
           if [ ! -f {GENOME}.fasta.fai ]; then
               samtools faidx {GENOME}.fasta
              fi
           """

      
rule gatk_recalibrate:
    input: JA_OUT+"bams/{x}.fin.bam", JA_OUT+"bams/{x}.fin.bam.bai"
    output: JA_OUT+"bams/{x}.fin.recal.bam"
    params: gatk=GATK,genome=GENOME+".fasta",rname="pl:hapcall",mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
#    threads: 1
	group: "align"
    shell: """
           {params.gatk} -T PrintReads -R {params.genome} -I {input[0]} -o {output[0]}  --defaultBaseQualities 40 --allow_potentially_misencoded_quality_scores
            """


rule picard_headers:
    input:  JA_OUT+"bams/{x}.sorted.bam"
    output: temp(JA_OUT+"bams/{x}.fin.bam")
    params: picard=PICARD,rname="pl:headers",mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100",rgpl=FOLDER
#    threads: 1
    group: "align"
    shell: """  
          {params.picard} I={input} O={output} RGID={wildcards.x} RGPL={params.rgpl} RGLB={wildcards.x} RGPU={wildcards.x} RGSM={wildcards.x} RGCN={wildcards.x} RGDS={input} Validation_Stringency=LENIENT
         """
         
         
rule fasta_dict:
    input: GENOME+".fasta"
    output: GENOME+".dict"
    params: picard=PICARD,rname="pl:headers",mem="16g",time="2:00:00",partition="ccr",gres="lscratch:100"
    shell: """
          module load picard;
          java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar CreateSequenceDictionary REFERENCE={input} OUTPUT={output}
           """

