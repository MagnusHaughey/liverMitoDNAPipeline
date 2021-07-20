
import numpy as np 
import sys

print(sys.argv)

outfile = sys.argv[1]
f = open(outfile , 'w')

# Write preamble
f.write('\\documentclass[10pt]{article}\n \\usepackage[english]{babel}\n \\usepackage{graphicx}\n \\usepackage[margin=1in]{geometry}\n \
	\\usepackage{xcolor, colortbl}\n \\usepackage{makecell}\n \\usepackage{pdfpages}\n \n\n\n \\begin{document}')

# Data name
f.write('\n\\section*{\\Large{' + outfile.split("/")[-2].replace('_' , '\\_') + '}}\n\n\n')


# Make table with flagstat and insert size metrics
# Flagstats
flagstats_file = open(sys.argv[2] , 'r')
f.write('\\vspace*{0.5in}\n \\footnotesize\n \\centering\n')
f.write('\\begin{tabular}{|c|c|}\n \\hline\n  \\cellcolor{gray!25} Samtools \\textit{flagstats} & \\makecell[l]{\\\\ \n')

for line in flagstats_file.readlines():
	f.write(line.replace('%' , '\\%').replace('>=' , '$\\geq$').replace('\n' , '\\\\ \n'))

f.write('\\vspace*{0.01in}} \\\\ \n\n')




# CollectInsertSizeMetrics
picard_inserts_file = open(sys.argv[3] , 'r')
f.write(' \\hline  \\cellcolor{gray!25} Picard \\textit{CollectInsertSizeMetrics} & \\makecell[l]{ \\\\ \n')

lines = np.genfromtxt(picard_inserts_file, skip_header=1, usecols=[i for i in range(18)], dtype=str)

for i in range(len(lines[0])):
	line = "{} = {}\n".format(lines[0][i] , lines[1][i])
	f.write(line.replace('_' , '\\_').replace('\n' , '\\\\ \n'))

f.write('\\vspace*{0.01in}} \\\\ \n\n')



# fastQC flags
fastqc_summary_file = open(sys.argv[4] , 'r')
f.write(' \\hline  \\cellcolor{gray!25} fastQC & \\makecell[l]{ \\\\ \n')

for line in fastqc_summary_file.readlines():
	if ('FAIL' in line):
		f.write('\\textbf{')
		f.write(line.replace('\n' , '}\\\\ \n'))
	else:
		f.write(line.replace('\n' , '\\\\ \n'))


f.write('\\vspace*{0.1in}} \\\\ \n \\hline\n \\end{tabular}\n\n \\clearpage \n\n')






# Then add figures and plots
for i in range(len(sys.argv)-3,len(sys.argv)):
#for i in range(6,9):
	path_fragments = sys.argv[i].split("/")
	relative_path = "./" + path_fragments[-2] + "/" + path_fragments[-1]
	f.write('\\clearpage \n \\includegraphics[width=1\\textwidth]{' + '{}'.format(relative_path) + '}\n\n')





f.write('\\end{document}')
f.close()



























