import os
import sys
import numpy as np
import time
import convert_ip

# some constants 
# CAUTION: DENSITY IS HARDWIRED AT THE MOMENT
BOLTZMANN = 1.38062e-16
EV2ERGS   = 1.602192e-12
LOG_N = 9 
COLUMN = 23
OUTPUT_FOLDER = "cld_output"


def run_cloudy(fname, cloudy_path="/Users/matthewsj/software/c17.01/source/"):
	isys = os.system("{}/cloudy.exe < {}.in > {}.out".format(cloudy_path, fname, fname))
	return (isys)

def get_param_string(root, log_xi, log_U, BB_temp = 580225.906078):
	'''
	Create a string to write to a cloudy input file. 

	Parameters:
		root 	str
				root filename to use for naming
		
		log_xi  float 
				log of xi ionization parameter (only used for title)

		log_U 	float 
				log of U ionization parameter
		
		BB_temp float 
				Blackbody temperature in K (linear)

	Returns:
		my_string str
		string containing text to use for input file
	'''
	my_string = '''title BB model, log density {:.1f}, log xi {:.1f}
# =========
# commands controlling continuum =========
Blackbody {:.1f}
ionization parameter {:.1f}
# =========
# commands for density & abundances
# =========
hden {:.1f}
# =========
# commands controlling geometry  
# =========
stop column density 19
# =========
# other commands for details
# =========
iterate convergence 
# these are to try to speed things up
init "c84.ini"
no molecules
no level2 lines
no fine opacities
atom h-like levels small
atom he-like levels small
COSMIC RAY BACKGROUND
element limit off -8
# =========
# commands controlling output    =========
# =========
normalize to "H  1" 1215.67A 100 
print line faint 1  # get rid of lots of faint lines c
# print off
print lines column
save lines, array "{}.lines" last
save continuum "{}.cont" 
#
# ========================================
'''.format(LOG_N, log_xi, BB_temp, log_U, LOG_N, root, root)
	return my_string

def initialise_cloudy_sim(root, log_xi, log_U, BB_temp = 580225.906078):
	'''
	Initialise a cloudy simulation 

	Parameters:
		root 	str
				root filename to use for naming
		
		log_xi  float 
				log of xi ionization parameter (only used for title)

		log_U 	float 
				log of U ionization parameter
		
		BB_temp float 
				Blackbody temperature in K (linear)

	Returns:
		my_string str
		string containing text to use for input file
	'''

	# get the string to write to file
	params = get_param_string(root, log_xi, log_U, BB_temp = BB_temp)
	
	# chuck the text in an input file 
	fname = "{}.in".format(root)
	f = open(fname, "w")
	f.write(params)
	f.close()

def run_logxi_grid(logxis =  np.arange(-3,2,0.5), BB_TEMP = 50.0):
	'''
	Run a grid of blackbody cloudy simulations 

	Parameters:	
		logxis	array-like
				array of floats containing log10 of ionization parameter xi
		
		BB_temp float 
				Blackbody temperature in ELECTRON VOLTS (linear)
	'''

	for i, logxi in enumerate(logxis):
		print (i, logxi)

		# we need to calculate U, the IP used by Cloudy. 
		xi = 10.0 ** logxi 
		factor = convert_ip.get_conversion_factor_blackbody(temp_eV = 50.0)
		U = xi * factor 
		logU = np.log10(U)

		# root filename to use for input and output files
		root = "BB_{:.1f}_1858_n{:.1f}_xi{:.1f}".format(BB_TEMP, LOG_N, logxi)

		print (BB_TEMP * EV2ERGS / BOLTZMANN)

		# Initialise simulation by creating input file 
		initialise_cloudy_sim(root, logxi, logU, BB_temp = BB_TEMP * EV2ERGS / BOLTZMANN)

		# run the sim! 
		run_cloudy(root)

if __name__ == "__main__":
	run_logxi_grid(logxis =  np.arange(-3,2,0.5), BB_TEMP = 50.0)

