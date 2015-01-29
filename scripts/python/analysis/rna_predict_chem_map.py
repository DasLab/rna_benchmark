#!/usr/bin/python

###############################################################################

import argparse
import numpy as np
import subprocess
from glob import glob 
import os
from os.path import basename, dirname, exists, isfile, isdir
import string
from sys import exit 

###############################################################################

parser = argparse.ArgumentParser(description='Predict DMS reactivity from pdb.')
#parser.add_argument('-target_files', nargs='+', help='List of additional target files.', default=['favorites.txt','favorites2.txt'])
parser.add_argument('-targets', nargs='+', help='List of targets.', default=['*'])
parser.add_argument('-silent_file', default='swm_rebuild.out')
parser.add_argument('-ndecoys', default=20)
parser.add_argument('-analysis_dir',default='analysis')
parser.add_argument('-path_to_rosetta_exe',default='')
args = parser.parse_args()

###############################################################################

'''
This script is a wrapper for rna_predict_chem_map, a Rosetta application

TODO:
1. cd < target directory > 
2. easy_cat.py SWM, if swm_rebuild.out not found, extract 20 lowest energy pdbs
3. Call rna_predict_chem_map on 20 lowest energy pdbs
4. Analyze output (.DMS files) and calculate statistics for target (swm_rebuild.out.10.pdb.DMS.txt)
5. Save statistics for target
6. cd ..
7. Plot prediction statistics/maps for all targets 

'''

###############################################################################

homedir = os.getcwd()
if args.targets == ['*']:
	targets = [dir for dir in glob('*') if isdir(dir)]
else:
	targets = args.targets  
silent_file = args.silent_file
ndecoys = int(args.ndecoys)
analysis_dir = args.analysis_dir
path_to_rosetta_exe = args.path_to_rosetta_exe
if len(path_to_rosetta_exe):
	assert(exists(path_to_rosetta_exe))
	path_to_rosetta_exe += '/'

###############################################################################

DMS_reactivity_per_target = {}
for target in targets:
	
	### print current target
	print '\n', target

	### go into target directory 
	os.chdir(target)
	print os.getcwd()

	### make analysis directory
	while not(exists(analysis_dir)):
		cmdline = ['mkdir','analysis']
		out,err = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if not(exists(analysis_dir)):
			print 'CMDLINE:', cmdline
			print 'OUT:', out
			print 'ERROR: ', err
	assert(exists(analysis_dir))
	os.chdir(analysis_dir)
	print os.getcwd()

	### cat silent_files
	while not(exists(silent_file)):
		cmdline = ['easy_cat.py','../SWM']
		out,err = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if not(exists(silent_file)):
			print 'CMDLINE:', cmdline
			print 'OUT:', out
			print 'ERROR:', err
	assert(exists(silent_file))

	### append target name to silent_file name
	target_silent_file = target + '_' + silent_file 
	while not(exists(target_silent_file)):
		cmdline = ['cp',silent_file,target_silent_file]
		out,err = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if not(exists(target_silent_file)):
			print 'CMDLINE:', cmdline
			print 'OUT:', out
			print 'ERROR:', err
	assert(exists(target_silent_file))

	### extract 20 lowest energy decoys
	target_pdbs = [target_silent_file + '.%d.pdb' % idx for idx in xrange(1,ndecoys+1)]
	assert(len(target_pdbs) == ndecoys)
	while not (sum([exists(f) for f in target_pdbs]) == len(target_pdbs)):
		cmdline = ['extract_lowscore_decoys.py', target_silent_file,str(ndecoys)] 
		out,err = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if not (sum([exists(f) for f in target_pdbs]) == len(target_pdbs)):
			print 'CMDLINE:', cmdline
			print 'OUT:', out
			print 'ERROR:', err
	assert(sum([exists(f) for f in target_pdbs]) == len(target_pdbs))
	print silent_file
	print target_silent_file
	print string.join(target_pdbs, '\n')


	### NOTE: might be smarter to wrap everything above in an init() function or separate 
	### post-processing script.
	### NOTE: could check for DMS files before running everything above

	### call rna_predict_chem_map on 20 lowest energy pdbs for target
	dms_txt_files = []
	for target_pdb in target_pdbs:
		dms_txt_file = target_pdb+'.DMS.txt'
		dms_txt_files.append(dms_txt_file)
		while not exists(dms_txt_file):
			cmdline = [path_to_rosetta_exe+'rna_predict_chem_map', '-s', target_pdb]
			out,err = subprocess.Popen(cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
			if not exists(dms_txt_file):
				print 'CMDLINE:', cmdline
				print 'OUT:', out
				print 'ERROR:', err
		assert(exists(dms_txt_file))
		print dms_txt_file	
	assert(sum([exists(f) for f in dms_txt_files]) == len(dms_txt_files)) 
	print string.join(dms_txt_files, '\n')


	### find average mean DMS reactivites at each residue
	DMS_mean_sum_per_res = {}
	DMS_mean_mean_per_res = {}
	for dms_txt_file in dms_txt_files:
		dms_txt_file_lines = open(dms_txt_file,  'r' ).readlines()
		col_names = dms_txt_file_lines[0].split()
		for idx,line in enumerate(dms_txt_file_lines[1:],start=1):
			cols = line.split()
			chain = cols[0]
			resid = cols[1]
			DMS_mean = float(cols[2])
			chain_resid = chain + ':' + resid
			if idx in DMS_mean_sum_per_res.keys():
				DMS_mean_sum_per_res[ idx ] += DMS_mean
			else:
				DMS_mean_sum_per_res[ idx ] = DMS_mean

	for idx, DMS_mean_sum in DMS_mean_sum_per_res.iteritems():
		DMS_mean_mean_per_res[ idx ] = DMS_mean_sum / len(dms_txt_files)

	DMS_reactivity_per_target[ target ] = DMS_mean_mean_per_res

	### return to run directory
	os.chdir( homedir )
	
for target, DMS_mean_mean_per_res in DMS_reactivity_per_target.iteritems():
	print '\n', target

	for idx, DMS_mean_per_res in DMS_mean_mean_per_res.iteritems():
		print str(idx)+':',DMS_mean_per_res

###
### PRINT HEAT MAP OF DATA
###

'''
(y = residues)
1  -|
2  -|
3  -|    *  *  *     *  *  *     *  *  *  *  *  *  *
4  -|
5  -| *  *     *  *  *  *  *  *  *  *  *  *  *  *  *
6  -| *  *  *  *  *  *     *  *  *  *     *  *     *
7  -|
8  -|
9  -|
10 -|
11 -| *     *  *  *  *  *  *  *  *  *     *  *  *  *
12 -|
13 -| *  *  *     *  *     *  *     *  *  *  *     *
14 -|
15 -|
16 -|________________________________________________
	  01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16
	  					(x = targets)

WHERE: * == mean represented as some color / value



'''















