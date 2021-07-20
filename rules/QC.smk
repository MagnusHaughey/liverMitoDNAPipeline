# An example collection of Snakemake rules imported in the main Snakefile.



def get_all_fastQC(wildcards):
        return sorted( glob.glob( "results/QC/fastQC/GC_EC_*" + wildcards.manifest + "_" + wildcards.patient + "_" + wildcards.sample + "_" + wildcards.primer + "_L*"))





rule perform_fastqc:
	input:
		bam = config['BAMfiles'] + "GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
	output:
		zip = "results/QC/fastQC/GC_EC_{manifest}_{patient}_{sample}_{primer}_fastqc.zip",
		html = "results/QC/fastQC/GC_EC_{manifest}_{patient}_{sample}_{primer}_fastqc.html",
	threads: 1
	shell:
		"""
		module load fastqc
		fastqc {input.bam} -o ./results/QC/fastQC/
		"""






rule summarise_fastqc:
	input:
		zip = "results/QC/fastQC/GC_EC_{manifest}_{patient}_{sample}_{primer}_fastqc.zip",
	output:
		summary = "results/QC/fastQC_summaries/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
	threads: 1
	shell:
		"""
		bash scripts/shell/getQC_fails.sh {input.zip} > {output.summary}
		"""







rule samtools_flagstat:
	input:
		bam = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam"
	output:
		flagstat = "results/QC/flagstats/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
	threads: 1
	shell:
		"""
		module load samtools
		samtools flagstat {input.bam} > {output.flagstat}
		"""







rule picard_qualityScore:
	input:
		bam = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
	output:
		txt = "results/QC/picardQualityScoreDistribution/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
		pdf = "results/QC/picardQualityScoreDistribution/GC_EC_{manifest}_{patient}_{sample}_{primer}.pdf"
	params:
		picard = config['picard'],
	shell:
		"""
		module load java
		module load R
		java -jar -Xmx4G {params.picard} QualityScoreDistribution \
 		 I={input.bam} \
 		 O={output.txt} \
 		 CHART={output.pdf}
		"""






rule picard_insertSize:
	input:
		bam = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
	output:
		txt = "results/QC/picardCollectInsertSizeMetrics/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
		pdf = "results/QC/picardCollectInsertSizeMetrics/GC_EC_{manifest}_{patient}_{sample}_{primer}.pdf"
	params:
		picard = config['picard'],
	shell:
		"""
		module load java
		module load R
		java -jar -Xmx4G {params.picard} CollectInsertSizeMetrics \
		 I={input.bam} \
		 O={output.txt} \
		 H={output.pdf}
		"""






rule combine_QC_summaries:
	input:
		fastQC_zip = "results/QC/fastQC/GC_EC_{manifest}_{patient}_{sample}_{primer}_fastqc.zip",
		fastQC_html = "results/QC/fastQC/GC_EC_{manifest}_{patient}_{sample}_{primer}_fastqc.html",
		fastQC_summary = "results/QC/fastQC_summaries/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
		flagstat = "results/QC/flagstats/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
		picardQScore_pdf = "results/QC/picardQualityScoreDistribution/GC_EC_{manifest}_{patient}_{sample}_{primer}.pdf",
		picardInsertSize_txt = "results/QC/picardCollectInsertSizeMetrics/GC_EC_{manifest}_{patient}_{sample}_{primer}.txt",
                picardInsertSize_pdf = "results/QC/picardCollectInsertSizeMetrics/GC_EC_{manifest}_{patient}_{sample}_{primer}.pdf",
		coverage = "results/coverage/samtools_depth/GC_EC_{manifest}_{patient}_{sample}_{primer}_cumulative_coverage.png",
		rawCoverage = "results/coverage/samtools_depth/GC_EC_{manifest}_{patient}_{sample}_{primer}_cumulative_coverage_rawCounts.png",
	output:
		summary_tex = "results/QC/summaries/GC_EC_{manifest}_{patient}_{sample}_{primer}/QC_summary.tex",
	shell:
		"""

			summary_dir="results/QC/summaries/GC_EC_"{wildcards.manifest}"_"{wildcards.patient}"_"{wildcards.sample}"_"{wildcards.primer}
			mkdir -p $summary_dir
	
			mkdir -p $summary_dir'/fastQC/' && cp {input.fastQC_zip} {input.fastQC_html} $summary_dir'/fastQC/'
			mkdir -p $summary_dir'/picardQualityScoreDistribution/' && cp {input.picardQScore_pdf} $summary_dir'/picardQualityScoreDistribution/'
			mkdir -p $summary_dir'/picardCollectInsertSizeMetrics/' && cp {input.picardInsertSize_pdf} $summary_dir'/picardCollectInsertSizeMetrics/'
			mkdir -p $summary_dir'/coverage/' && cp {input.coverage} $summary_dir'/coverage/'
			cp {input.rawCoverage} $summary_dir'/coverage/'

			grep "## METRICS" -A 2 {input.picardInsertSize_txt} >> {input.picardInsertSize_txt}.tmp

			picardQS_chart=$summary_dir"/picardQualityScoreDistribution/GC_EC_"{wildcards.manifest}"_"{wildcards.patient}"_"{wildcards.sample}"_"{wildcards.primer}".pdf"
			picardInsertSize_chart=$summary_dir"/picardCollectInsertSizeMetrics/GC_EC_"{wildcards.manifest}"_"{wildcards.patient}"_"{wildcards.sample}"_"{wildcards.primer}".pdf"
			coverage=$summary_dir"/coverage/GC_EC_"{wildcards.manifest}"_"{wildcards.patient}"_"{wildcards.sample}"_"{wildcards.primer}"_cumulative_coverage.png"
			rawCoverage=$summary_dir"/coverage/GC_EC_"{wildcards.manifest}"_"{wildcards.patient}"_"{wildcards.sample}"_"{wildcards.primer}"_cumulative_coverage_rawCounts.png"


			python3 scripts/python/make_QCsummaries.py {output.summary_tex} {input.flagstat} {input.picardInsertSize_txt}.tmp {input.fastQC_summary} $picardQS_chart $picardInsertSize_chart $coverage $rawCoverage

		"""






