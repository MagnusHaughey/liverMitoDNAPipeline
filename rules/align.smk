

rule sort_bam:
        input:
                bam = config['BAMfiles'] + "GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
        output:
                bam = "results/BAM/2.sorted/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
        params:
                picard = config["picard"]
        threads: 1
	shell:
                """
                module load java
                java -jar -Xmx4G {params.picard} SortSam \
                 INPUT={input.bam} \
                 OUTPUT={output.bam} \
                 SORT_ORDER=coordinate
                """


rule markDuplicates_index:
	input:
		bam = "results/BAM/2.sorted/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
	output:
		bam = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
		bai = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bai",
		metrics = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}_MarkDuplicatesMetrics.txt",
	params:	
		picard = config['picard'],
	threads: 1
	shell:
		"""
		module load java
        	java -jar -Xmx4G {params.picard} MarkDuplicates \
            	 I={input.bam} \
            	 O={output.bam} \
            	 M={output.metrics} \
            	 CREATE_INDEX=true \
		 REMOVE_DUPLICATES=true
		"""

