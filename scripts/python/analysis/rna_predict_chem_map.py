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
reactive_residues_per_target = {}
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
	reactive_residues = {}
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

			if idx in reactive_residues.keys():
				reactive_residues[ idx ] += ( DMS_mean > 1.0 ) 
			else:
				reactive_residues[ idx ] = ( DMS_mean > 1.0 )
			
	for idx, DMS_mean_sum in DMS_mean_sum_per_res.iteritems():
		DMS_mean_mean_per_res[ idx ] = DMS_mean_sum / len(dms_txt_files)

	DMS_reactivity_per_target[ target ] = DMS_mean_mean_per_res
	reactive_residues_per_target[ target ] = reactive_residues




	### return to run directory
	os.chdir( homedir )

for target, DMS_mean_mean_per_res in DMS_reactivity_per_target.iteritems():
	print '\n', target
	for idx, DMS_mean_per_res in DMS_mean_mean_per_res.iteritems():
		print str(idx)+':',DMS_mean_per_res

###############################################################################

data = []
residues = sorted(DMS_reactivity_per_target[ targets[0] ].keys())

for res_idx in residues:
	data_row = []
	for target in sorted(DMS_reactivity_per_target):
		data_row.append(DMS_reactivity_per_target[ target ][ res_idx ])
	data.append(np.array(data_row))
data = np.array(list(reversed(data)))

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

rows = list(reversed(residues))
columns = [ x for x in xrange(1, len(DMS_reactivity_per_target.keys())+1)]

print 'DATA: ', data
print 'ROWS: ', rows
print 'COLS: ', columns

save_fig_dir = '../../../Figures/'
run_dir = basename(os.getcwd())
fullpdfname = save_fig_dir + 'DMS_reactivity_predictions_'+run_dir+('_%d_decoys'%ndecoys)+'.pdf'

print '\nMaking figure in: %s\n' % fullpdfname
pp = PdfPages( fullpdfname )


fig,ax = plt.subplots()
fig.set_size_inches(11, 8.5)

plt.pcolor(data, cmap=plt.get_cmap('Blues'))
plt.colorbar( orientation='horizontal')
plt.xticks(np.arange(0,len(columns))+0.5)
plt.yticks(np.arange(0,len(rows))+0.5)

ax.set_xlim(0,len(columns))
ax.set_ylim(0,len(rows))

ax.xaxis.tick_top()
ax.yaxis.tick_left()

ax.set_xticklabels(columns, minor=False, fontsize=12)
ax.set_yticklabels(rows, minor=False, fontsize=12)
ax.set_aspect( 1 ) #'equal')

plt.text(0.5,1.25, 'DMS Reactivity Predictions\n\n%s' % run_dir,
		 fontsize=14,
		 horizontalalignment='center',
		 transform=ax.transAxes
	    )
plt.ylabel('Residue', fontsize=12)
plt.xlabel('Target', fontsize=12)


# save as pdf and close
pp.savefig()
pp.close()

# open pdf
out, err = subprocess.Popen(['uname'], stdout=subprocess.PIPE).communicate()
if 'Darwin' in out:
	subprocess.call(['open',fullpdfname])
if 'Linux' in out:
	subprocess.call(['xdg-open',fullpdfname])


for target in sorted(reactive_residues_per_target):
	print '\n', target
	reactive_residues = reactive_residues_per_target[ target ]

	for idx, reactivity in reactive_residues.iteritems():
		if reactivity > 0:
			print '%d: %d/%d' % (idx,reactivity,ndecoys)