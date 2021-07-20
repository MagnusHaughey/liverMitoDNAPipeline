


rule samtools_depth:
	input:
		bam = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
	output:
		depth_plot = "results/coverage/samtools_depth/GC_EC_{manifest}_{patient}_{sample}_{primer}_cumulative_coverage.png",
		raw_coverage_plot = "results/coverage/samtools_depth/GC_EC_{manifest}_{patient}_{sample}_{primer}_cumulative_coverage_rawCounts.png",
	threads: 1
	shell:
		"""
		module load samtools
		
		samtools depth -q 20 -d 25000 {input.bam} | python3 scripts/python/plot_cumulative_coverage.py \
										-O {output.depth_plot} \
										--dataType RAW_COUNT \
										--primer {wildcards.primer} >> results/coverage/all_depths.txt
		"""


