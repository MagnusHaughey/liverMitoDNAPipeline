
import sys

# Ensure python can see numpy etc. packages
#sys.path.append("/data/BCI-EvoCa2/magnus/liver_mtDNA/mtDNApipeline/.venv/lib/python3.6/site-packages/")


import numpy as np
import matplotlib.pyplot as plt
from pylab import text
import argparse
from operator import itemgetter
import math



# Parse command line arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-I', help='Input file with raw coverage data')
parser.add_argument('-O', help='File path for outputted .png', type=str)
parser.add_argument('--dataType', help='Specify if data is in format of "coordinate, depth" or "bin, count"', type=str)
parser.add_argument('--primer', help='Specify if primer M1 or M2, default is both M1+M2 combined', type=str, default='')
args = parser.parse_args()


if not ( (args.dataType == 'BINNED') or (args.dataType == 'RAW_COUNT') ): 
	raise Exception('--dataType argument should be either "BINNED" or "RAW_COUNT".')

if not ( (args.primer == 'M1A') or (args.primer == 'M1B') or (args.primer == 'M2A') or (args.primer == 'M2B') or (args.primer == '') ): 
	raise Exception('--primer argument should be either "M1*", "M2*", or not given.')


# Obtain patient, sample and primer variables from output 
params = args.O.split("/")[-1]
params = params.strip("GC_EC_")
params = params.strip("_cumulative_coverage.png")
manifest, patient, sample, primer = params.split("_")



# Import data from stdin
coverage_raw_data = sys.stdin.read()
coverage_raw_data = coverage_raw_data.splitlines()



# Read depth at each position from third column data into ' depth[position , depth] '
position = []
depth = []
bins = []


if (args.dataType == 'RAW_COUNT'):
	for i in range(len(coverage_raw_data)):
			line_split = coverage_raw_data[i].split("\t")
			position.append(int(line_split[1]))	
			depth.append(int(line_split[2]))

else:
	for i in range(1,len(coverage_raw_data)-1):
			line_split = coverage_raw_data[i].split("\t")
			bins.append(int(line_split[0]))
			depth.append(int(line_split[1]))




if (args.dataType == 'RAW_COUNT'):

	# Find approximate on/off amplicon coordinates
	switch_position_1 = 0
	switch_position_2 = 0
	primer_on = (depth[0] >= 100)
	for i in range(len(depth)):
		if (depth[i] < 100) and (primer_on == True):
			switch_position_1 = position[i]                # primer off position
			primer_on = False

		if (depth[i] > 100) and (primer_on == False):
			switch_position_2 = position[i]                # primer on position
			primer_on = True

		if (switch_position_1 != 0) and (switch_position_2 != 0):
			break


	plt.switch_backend('agg')

	plt.rc('xtick',labelsize=14)
	plt.rc('ytick',labelsize=14)
	plt.xticks([0,switch_position_1,switch_position_2])

	plt.plot(position , depth, color='blue', linewidth=0.75)
	axes = plt.gca()
	axes.set_xlim([0,16571])
	plt.xlabel('Position', fontsize=16)
	plt.ylabel('#reads', fontsize=16)

	# Draw vertical lines at approximate on/off coordinates
	#plt.axvline(x=switch_position_1)
	#plt.axvline(x=switch_position_2)

	plt.tight_layout()


	# Export plot
	out_fig = args.O.strip(".png") + "_rawCounts.png"
	plt.savefig(out_fig , dpi=300 , format='png')





# Generate binned data 
if (args.dataType == 'RAW_COUNT'):
	bins = np.linspace(0 , 25000 , 25000)
	
	# if only one (i.e. either M1 or M2) specified, then only use depth data for the appropriate region
	if ('M1' in args.primer):
		depth = [val for val in depth if (val > 500)]
		mean_depth = np.mean(depth)
	if ('M2' in args.primer):
		depth = [val for val in depth if (val > 500)]
		mean_depth = np.mean(depth)

	depth, bins = np.histogram(depth , bins=bins)
	bins = bins[:-1]


else: 
	# Picard counts coverage for whole genome. So cut first ~10 bins to avoid all the zero counts and counts from mis-aligned reads
	for i in range(500):
		depth[i] = 0
		
	depth = np.asarray(depth)


# Occasionally sequencing fails and bam file has no high quality reads. In this case, print nothing and output empty graph 
# Otherwise remove trailing zeros from depth list, and corresponding bin value
while (len(depth) > 0) and (depth[-1] == 0):
	depth = depth[:-1]
	bins = bins[:-1]




cumulative_coverage = np.zeros(len(bins))
if (args.dataType == 'BINNED'):
	mean_depth = 0.0

for i in range(len(bins)):
	cumulative_coverage[i] = np.sum(depth[i:])
	if (args.dataType == 'BINNED'):
		mean_depth += depth[i]*bins[i]

if (args.dataType == 'BINNED'):
	mean_depth /= np.sum(depth)




# Normlise
if (len(cumulative_coverage) > 0):
	cumulative_coverage /= (np.max(cumulative_coverage)/100.0)


#============= Plotting stuff

plt.switch_backend('agg')

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.xticks([0,5000,10000,15000,20000,25000])

plt.plot(bins , cumulative_coverage, color='blue')
axes = plt.gca()
axes.set_xlim([0,25000])
axes.set_ylim([0,100])
plt.xlabel('Depth of coverage', fontsize=16)
plt.ylabel('% of target', fontsize=16)
plt.title('Cumulative coverage', fontsize=16)

if (math.isnan(mean_depth) == False):
	text(0.70, 0.90,'Mean depth = {}'.format(int(mean_depth)), fontsize=16, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes)
	print("{} {} {} {}".format(patient, sample, primer, int(mean_depth)))
else:
	text(0.70, 0.90,'Mean depth = {}'.format("-"), fontsize=16, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes)
	print("{} {} {} {}".format(patient, sample, primer, '-'))

plt.fill_between(bins , cumulative_coverage, y2=0)
plt.tight_layout()

# Export plot
out_fig = args.O
plt.savefig(out_fig , dpi=300 , format='png')



# Finally, print the mean depth to stdout 
#print("{} {} {} {}".format(patient, sample, primer, int(mean_depth)))
#print("{}".format(int(mean_depth)))

















