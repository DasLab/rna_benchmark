#!/usr/local/bin/python

###############################################################################

import argparse
import numpy as np
import subprocess
from glob import glob
import os
from os.path import basename, dirname, exists, isfile, isdir
import string
from sys import exit
from rdat import *

###############################################################################

parser = argparse.ArgumentParser(description='Process DMS chemical mapping data.')
parser.add_argument('-exp_rdat_files', nargs='*', help='RDAT files for processing', default=[])
parser.add_argument('-pred_rdat_files', nargs='*', help='RDAT files for processing', default=[])
parser.add_argument('-exp_seqpos', nargs='*', help='list of sequence positions to compare', default=[])
parser.add_argument('-pred_seqpos', nargs='*', help='list of sequence positions to compare', default=[])

args = parser.parse_args()

###############################################################################

rdat = RDATFile()

if len(args.exp_rdat_files) and len(args.pred_rdat_files):

	experimental_data = {}
	prediction_data = {}

	for rdat_file in args.exp_rdat_files:
		experimental_data[ rdat_file ] = rdat.get_reactivity(rdat_file)

	for rdat_file in args.pred_rdat_files:
		prediction_data[ rdat_file ] = rdat.get_reactivity(rdat_file)

	pdfname = 'experimental_vs_prediction_scatter_plot.pdf'
	title = basename( pdfname.replace('_',' ').replace('.pdf','') )

	scatter_plot_experimental_vs_prediction( experimental_data, prediction_data, args.exp_seqpos, args.pred_seqpos, fullpdfname=pdfname, title=title )