
import numpy as np 
import sys
import argparse


# Parse command line arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-I', help='Input file with raw coverage data')
parser.add_argument('--input_one', help='', type=str)
parser.add_argument('--input_two', help='', type=str)
parser.add_argument('--input_three', help='', type=str)
parser.add_argument('--output', help='', type=str)
args = parser.parse_args()


# Read in data files
position , p_val , raw_freq = np.loadtxt(args.input_one , unpack=True , usecols=(2,5,6))
ref_base , var_base = [] , []

for line in open(args.input_one , 'r').readlines():
	fields = line.replace('   ',' ').replace('  ' , ' ' ).split(" ")
	ref_base.append(fields[3])
	var_base.append(fields[4])



n_tst_fw , cov_tst_fw , n_tst_bw , cov_tst_bw , n_ctrl_fw = np.loadtxt(args.input_two , unpack=True , skiprows=1 , usecols=(2,3,4,5,6))
cov_ctrl_fw , n_ctrl_bw , cov_ctrl_bw = np.loadtxt(args.input_three , unpack=True , skiprows=1 , usecols=(1,2,3))


# Compute "shifted" variant frequencies
shifted_var_freq = (( n_tst_fw + n_tst_bw )/( cov_tst_fw + cov_tst_bw ))



# Filtering 
filtered_out = []
f = open(args.output + '.METRICS.dat' , 'w')

if not isinstance(n_tst_bw, np.float64):
	for i in range(len(n_tst_bw)):

		# Filter on raw number of calls for each variant
		if ((n_tst_bw[i] + n_tst_fw[i]) < 10):
			filtered_out.append(i)
			if ((n_ctrl_bw[i] + n_ctrl_fw[i]) <= 10):
				f.write("Removed somatic mutation {}{}{} due to small number of raw calls\n".format(int(position[i]) , ref_base[i] , var_base[i]))
			else:
				f.write("Removed germline mutation {}{}{} due to small number of raw calls\n".format(int(position[i]) , ref_base[i] , var_base[i]))


elif isinstance(n_tst_bw, np.float64):
	
	# Filter on raw number of calls for each variant
	if ((n_tst_bw + n_tst_fw) < 10):
		filtered_out.append(0)
		if ((n_ctrl_bw + n_ctrl_fw) <= 10):
			f.write("Removed entry {}{}{} due to small number of raw calls\n".format(int(position) , ref_base , var_base))
		else:
			f.write("Removed germline mutation {}{}{} due to small number of raw calls\n".format(int(position) , ref_base , var_base))

f.close()




# Open output file 
if not isinstance(n_tst_bw, np.float64):

	#g = open(args.output , 'w')
	#
	#for i in range(len(shifted_var_freq)):
	#	if (i in filtered_out):
	#		continue
	#	else:
	#		g.write("{}{}{} {:1.10f} {}\n".format(int(position[i]) , ref_base[i] , var_base[i] , shifted_var_freq[i] , p_val[i]))
	#
	#g.close()


	# Write somatic calls to file
	g = open(args.output , 'w')
	for i in range(len(shifted_var_freq)):


		# If variant detected in control at frequency greater than 1%, define as germline
		if (( n_ctrl_fw[i] + n_ctrl_bw[i] )/( cov_ctrl_fw[i] + cov_ctrl_bw[i] ) <= 0.01):
		 	g.write("{}{}{} {:1.10f} {}\n".format(int(position[i]) , ref_base[i] , var_base[i] , shifted_var_freq[i] , p_val[i]))


	g.close()



else:

	#g = open(args.output , 'w')
	#
	#if (len(filtered_out) == 0):
	#	g.write("{}{}{} {:1.10f} {}\n".format(int(position) , ref_base , var_base , shifted_var_freq , p_val))
	#
	#
	#g.close()


	# Write somatic calls to file
	g = open(args.output , 'w')

	# If variant detected in control at frequency greater than 1%, define as germline
	if (len(filtered_out) == 0) and (( n_ctrl_fw + n_ctrl_bw )/( cov_ctrl_fw + cov_ctrl_bw ) <= 0.01):
		g.write("{}{}{} {:1.10f} {}\n".format(int(position) , ref_base , var_base , shifted_var_freq , p_val))

	g.close()














