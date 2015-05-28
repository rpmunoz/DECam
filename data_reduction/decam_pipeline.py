#! /usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import sys,os
import os.path
import time
import subprocess
import numpy as np
import pyfits
import multiprocessing, Queue
import ctypes
import matplotlib.pyplot as plt
import scipy.interpolate
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils import detect_sources, segment_properties, properties_table
from photutils.background import Background
from scipy.interpolate import Rbf

from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)

def main(argv):
	recipe =argv[0]
	input_file = ''
	output_file = ''
	try:
		opts, args = getopt.getopt(argv[1:],"hi:o:",["input=","output="])
	except getopt.GetoptError:
		print 'decam_pipeline.py -i <input_file> -o <output_file>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'decam_pipeline.py -i <input_file> -o <output_file>'
			sys.exit()
		elif opt in ("-i", "--input"):
			input_file = arg
		elif opt in ("-o", "--output"):
			output_file = arg
	print 'Recipe is "', recipe
	print 'Input file is "', input_file
	print 'Output file is "', output_file

if __name__ == "__main__":
	main(sys.argv[1:])

	n_cpu=2
	n_core=6
	n_processes=n_cpu*n_core*1


