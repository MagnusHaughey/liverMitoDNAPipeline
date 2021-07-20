
suppressMessages(library(deepSNV))
suppressMessages(library(argparse))


# Parse command line args
parser <- ArgumentParser()
#parser$add_argument('--binsize', type = 'double',
 #                   help="Binsize for deepSNV region specifying", default = 500)
parser$add_argument('--test', type = 'character',
                    help="Specify input secuencing data in .bam format")
parser$add_argument('--control', type = 'character',
                    help="Specify reference genome")
parser$add_argument('--primer', type = 'character',
		    help="Specify primer M1 or M2, or leave blank for the whole of the mtDNA")
parser$add_argument('--outpath', type = 'character',
                    help="Specify path for outputted PDF scatter plot")
args <- parser$parse_args()



regions <- data.frame(chr="gi|251831106|ref|NC_012920.1|", start = 1, stop=16571)


variants <- deepSNV(test = args$test, control = args$control, regions=regions, q=20)


outfile <- file.path(paste(args$outpath, sep = ""))
pdf(file = outfile)
	plot(variants)	# scatter plots
dev.off()



#show(variants) # show method
#summary(variants) # summary with significant SNVs
SNVs <- summary(variants, sig.level=0.05, adjust.method="BH")
show(SNVs)










