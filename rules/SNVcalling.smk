

rule deepSNV:   
	input:
		bam = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_{sample}_{primer}.bam",
		control = "results/BAM/3.indexed/GC_EC_{manifest}_{patient}_Stroma_{primer}.bam",
	output:
		scatterPlot = "results/SNVs/GC_EC_{manifest}_{patient}_{sample}_{primer}_scatterPlot.pdf",
		rawData = "results/SNVs/GC_EC_{manifest}_{patient}_{sample}_{primer}_rawDeepSNV.dat",
	shell:
		"""
		module load R/3.5.3

		if [ "{wildcards.sample}" != "Stroma" ];
		then

			Rscript scripts/R/execDeepSNV.R --test {input.bam} \
							--control {input.control} \
							--outpath {output.scatterPlot} \
							--primer {wildcards.primer} > {output.rawData}
		
		else
			touch {output.scatterPlot} {output.rawData}
		fi


		"""








rule summarise_SNV_reference_rCRS:
        input:
                rawData = "results/SNVs/GC_EC_{manifest}_{patient}_{sample}_{primer}_rawDeepSNV.dat",
        output:
                SNV_summary = "results/SNVs/GC_EC_{manifest}_{patient}_{sample}_{primer}_summary.dat",
        shell:
                """

                        size=$(stat --printf="%s" {input.rawData})
                        if (( $size > 342 ));
                        then

                                num_sig_variants_rep_A=$(grep "gi|251831106|ref|NC_012920.1|"  {input.rawData} | wc -l)
                                grep "gi|251831106|ref|NC_012920.1|"  {input.rawData} > {input.rawData}.temp_one
                                grep "n.tst.fw" -A $num_sig_variants_rep_A {input.rawData} > {input.rawData}.temp_two
				grep "cov.ctrl.fw" -A $num_sig_variants_rep_A {input.rawData} > {input.rawData}.temp_three

                                python3 scripts/python/summariseSNVs_rCRS.py --input_one {input.rawData}.temp_one --input_two {input.rawData}.temp_two --input_three {input.rawData}.temp_three --output {output.SNV_summary} 

                                rm {input.rawData}.temp_one {input.rawData}.temp_two {input.rawData}.temp_three

                        else

                                touch {output.SNV_summary}

                        fi


                """


