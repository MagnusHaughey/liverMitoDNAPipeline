# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

### Process config file
import glob
import os

configfile: "config.yaml"


#====================================


#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"


(MANIFEST, PATIENTNAMES, SAMPLENAMES, PRIMERS) = glob_wildcards( config['BAMfiles'] + "GC_EC_{manifest}_{patient}_{sample}_{primer}.bam" )



TARGETS = ["{manifest}_{patient}_{sample}_{primer}".format(manifest=manifest, patient=patient, sample=sample, primer=primer) for manifest, patient, sample, primer in zip(MANIFEST, PATIENTNAMES, SAMPLENAMES, PRIMERS)]
TARGETS = list(set(TARGETS))


#====================================


rule all:
	input:
		expand("results/QC/summaries/GC_EC_{target}/QC_summary.tex", target = TARGETS),
		expand("results/SNVs/GC_EC_{target}_summary.dat", target = TARGETS),



#====================================


include: "rules/QC.smk"
include: "rules/align.smk"
include: "rules/depth.smk"
include: "rules/SNVcalling.smk"
