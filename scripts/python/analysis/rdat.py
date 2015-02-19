#!/usr/local/bin/python

###############################################################################

import os
import subprocess
from os.path import basename, dirname, exists, isfile, isdir
import numpy as np
import string
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


###############################################################################

class RDATFile:
	def __init__(self):
		self.filename = None
		self.lines = None
		self.version = None
		self.name = None
		self.sequence = None
		self.structure = None
		self.annotation = None
		self.seqpos = {}
		self.reactivity = {}
		self.reactivity_error = {}
		self.header = []

	def clear(self):
		self.__init__()

	def load(self, filename):
		self.clear()
		self.filename = filename
		file = open( self.filename, 'r' )
		for line in file.readlines():
			line = line.replace('\n','')
			if not line:  
				continue
			if 'RDAT_VERSION' in line:
				self.version = line.split()[1:]
				self.header.append(line)
				continue
			if 'NAME' in line:
				self.name = line.split()[1:]
				self.header.append(line)
				continue
			if 'SEQUENCE' in line:
				self.sequence = line.split()[1:]
				self.header.append(line)
				continue
			if 'STRUCTURE' in line:
				self.structure = line.split()[1:]
				self.header.append(line)
				continue
			if 'ANNOTATION' in line:
				#self.annotation = line.split()[1:]
				self.header.append(line)
				continue
			if 'SEQPOS' in line:
				seqpos_list = line.split()[1:]
				for idx, seqpos in enumerate(seqpos_list, start=1):
					self.seqpos[ idx ] = seqpos.lower()
				self.header.append(line)
				continue
			if 'REACTIVITY:' in line:
				sample_idx = int(line.split()[0].split(':')[-1])
				reactivity_list = line.split()[1:]
				self.reactivity[ sample_idx ] = {}
				for idx, reactivity in enumerate(reactivity_list, start=1):
					self.reactivity[ sample_idx ][ idx ] = float(reactivity)
				continue
			if 'REACTIVITY_ERROR:' in line:
				sample_idx = int(line.split()[0].split(':')[-1])
				reactivity_error_list = line.split()[1:]
				self.reactivity_error[ sample_idx ] = {}
				for idx, reactivity_error in enumerate(reactivity_error_list, start=1):
					self.reactivity_error[ sample_idx ][ idx ] = float(reactivity_error)
				continue

	def write_data(self, filename, data):
		fout = open(filename, 'w')
		fout.write(string.join(self.header,'\n') + '\n\n')
		for sample_idx in sorted(data.keys()):
			fout.write('REACTIVITY:%s' % str(sample_idx) )
			for seqpos_idx in sorted( data[ sample_idx ].keys()):
				fout.write( '\t%f' % data[ sample_idx ][ seqpos_idx ] )
			fout.write('\n')
		fout.close()


	def print_(self):
		print self.filename
		print self.version
		print self.name
		print self.sequence
		print self.structure
		print self.annotation
		print self.seqpos
		print self.reactivity
		print self.reactivity_error


	def total_seqpos_reactivity(self, rdat_files):
		reactivity_totals = {}
		for rdat_file in rdat_files:
			reactivity = self.get_reactivity(rdat_file)
			for sample_idx in reactivity.keys():
				if not sample_idx in reactivity_totals.keys():
					reactivity_totals[ sample_idx ] = {}
				for seqpos_idx in reactivity[ sample_idx ].keys():
					if not seqpos_idx in reactivity_totals[ sample_idx ].keys():
						reactivity_totals[ sample_idx ][ seqpos_idx ] = 0.0
					reactivity_totals[ sample_idx ][ seqpos_idx ] += reactivity[ sample_idx ][ seqpos_idx ]
		return reactivity_totals

	def mean_seqpos_reactivity(self, rdat_files, return_sample_idx=None, verbose=False):
		mean_seqpos_reactivity = self.total_seqpos_reactivity(rdat_files)					
		for sample_idx in sorted(mean_seqpos_reactivity.keys()):
			for seqpos_idx in sorted(mean_seqpos_reactivity[ sample_idx ].keys()):
				mean_seqpos_reactivity[ sample_idx ][ seqpos_idx ] = mean_seqpos_reactivity[ sample_idx ][ seqpos_idx ] / len(rdat_files)
			if verbose:
				print 'SAMPLE:', sample_idx
				print 'NDECOYS:', len(rdat_files)
				print 'SEQPOS_IDX  MEAN_REACTIVITY'
				for seqpos_idx in list(reversed(sorted(mean_seqpos_reactivity[ sample_idx ].keys()))):
					print '%-10d  %f' %( seqpos_idx, mean_seqpos_reactivity[ sample_idx ][ seqpos_idx ] )
		if return_sample_idx:
			return mean_seqpos_reactivity[ return_sample_idx ]
		return mean_seqpos_reactivity

	def get_reactivity(self, rdat_file, seqpos_list=None):
		self.load(rdat_file)
		reactivity = self.reactivity
		if seqpos_list:
			for sample_idx in reactivity.keys():
				for seqpos_idx in reactivity[ sample_idx ].keys():
					if seqpos_idx in seqpos_list: continue
					del( reactivity[ sample_idx ][ seqpos_idx ] )
		return reactivity

	def get_reactivity_error(self, rdat_file, seqpos_list=None):
		self.load(rdat_file)
		reactivity_error = self.reactivity_error
		if seqpos_list:
			for sample_idx in reactivity_error.keys():
				for seqpos_idx in reactivity_error[ sample_idx ].keys():
					if seqpos_idx in seqpos_list: continue
					del( reactivity_error[ sample_idx ][ seqpos_idx ] )
		return reactivity_error

###############################################################################
# PLOTTING UTIL
###############################################################################

def open_figure(fullpdfname):
	out, err = subprocess.Popen(['uname'], stdout=subprocess.PIPE).communicate()
	if 'Darwin' in out:
		subprocess.call(['open',fullpdfname])
	if 'Linux' in out:
		subprocess.call(['xdg-open',fullpdfname])

def make_rna_chemical_map(reactivity_data, fullpdfname=None, title=None, ndecoys=None):

	sample_idx_list = sorted(reactivity_data.keys())
	seqpos_idx_list = sorted(reactivity_data[ sample_idx_list[0] ].keys())	
	
	data = []
	for seqpos_idx in seqpos_idx_list:
		seqpos_data = [ reactivity_data[ x ][ seqpos_idx ] for x in sample_idx_list ]
		data.append(np.array(seqpos_data)) 
	data = np.array(data)

	rows = seqpos_idx_list
	columns = [x for x in xrange(1, len(reactivity_data.keys())+1)]

	#print 'DATA: ', data
	#print 'ROWS: ', rows
	#print 'COLS: ', columns

	if not fullpdfname:
		save_fig_dir = './'#'../../../Figures/'
		run_dir = basename(os.getcwd())
		fullpdfname = '%s/DMS_reactivity_predictions_%s.pdf' % (save_fig_dir, run_dir)
		if ndecoys:
			fullpdfname.replace('.pdf', '_%d_decoys.pdf'%ndecoys)
	print '\nMaking figure in: %s\n' % fullpdfname
	
	pp = PdfPages( fullpdfname )

	fig,ax = plt.subplots()
	fig.set_size_inches(11, 8.5)

	plt.pcolor(data, cmap=plt.get_cmap('Greys'))
	plt.clim(0,2)
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

	if not title:
		title = 'DMS Reactivity Predictions\n\n%s' % run_dir
	plt.text(0.5,1.25, title,
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
	open_figure(fullpdfname)

def plot_seqpos_error(reactivity_data, seqpos1, seqpos2, fullpdfname=None, title=None, ndecoys=None):

	sample_idx_list = sorted(reactivity_data.keys())
	seqpos_idx_list = sorted(reactivity_data[ sample_idx_list[0] ].keys())	

	seqpos1_data = []
	seqpos2_data = []
	for sample_idx in sample_idx_list:
		seqpos1_data.append(reactivity_data[ sample_idx ][ seqpos1 ])
		seqpos2_data.append(reactivity_data[ sample_idx ][ seqpos2 ])

	seqpos1_data = np.array(seqpos1_data)
	seqpos2_data = np.array(seqpos2_data)

	#print 'seqpos1_data: ', seqpos1_data
	#print 'seqpos2_data: ', seqpos2_data

	if not fullpdfname:
		save_fig_dir = './' #'../../../Figures/'
		run_dir = basename(os.getcwd())
		fullpdfname = '%s/Error_Plot_Seqpos_%s_vs_Seqpos_%s_%s.pdf' % (save_fig_dir, str(seqpos2), str(seqpos1), run_dir)
		if ndecoys:
			fullpdfname.replace('.pdf', '_%d_decoys.pdf'%ndecoys)
	print '\nMaking figure in: %s\n' % fullpdfname
	
	pp = PdfPages( fullpdfname )

	fig,ax = plt.subplots()
	fig.set_size_inches(11, 8.5)

	plt.plot( seqpos1_data, seqpos2_data, marker='.', markersize=4, linestyle=' ' )
	plt.xlim(0,3)
	plt.ylim(0,3)
	plt.plot( plt.xlim(), plt.ylim(), color='black')

	if not title:
		title = 'Seqpos %s vs Seqpos %s Error\n\n%s' % (str(seqpos2), str(seqpos1), run_dir)
	plt.text(0.5,1.05,title,
			 fontsize=14,
			 horizontalalignment='center',
			 transform=ax.transAxes
		    )

	plt.ylabel('Seqpos %s' % str(seqpos2), fontsize=12)
	plt.xlabel('Seqpos %s' % str(seqpos1), fontsize=12)

	# save as pdf and close
	pp.savefig()
	pp.close()

	# open pdf
	open_figure(fullpdfname)


def scatter_plot_experimental_vs_prediction( experimental_data, prediction_data, exp_seqpos_list, pred_seqpos_list, fullpdfname=None, title=None ):

	print '\nMaking figure in: %s\n' % fullpdfname
	pp = PdfPages( fullpdfname )

	fig = plt.figure()

	# one plot per seqpos
	for idx, pred_seqpos in enumerate(pred_seqpos_list):

		exp_seqpos = exp_seqpos_list[ idx ]

		ax = fig.add_subplot( len(pred_seqpos_list), 1, idx+1 )

		print 'pred_seqpos =', pred_seqpos
		print 'exp_seqpos = ', exp_seqpos 

		for pred_fname in prediction_data.keys():
			pred_data = prediction_data[ pred_fname ]

			print 'pred_fname = ', pred_fname

			pred_sample_idx_list = sorted(pred_data.keys())

			pred_seqpos_data = []
			for sample_idx in pred_sample_idx_list:
				pred_seqpos_data.append(pred_data[ sample_idx ][ pred_seqpos ])
			pred_seqpos_data = np.array(pred_seqpos_data)


			for exp_fname in experimental_data.keys():
				exp_data = experimental_data[ exp_fname ]

				print 'exp_fname = ', exp_fname
				
				exp_sample_idx_list = sorted(exp_data.keys())

				exp_seqpos_data = []
				for sample_idx in exp_sample_idx_list:
					exp_seqpos_data.append(exp_data[ sample_idx ][ exp_seqpos ])
				exp_seqpos_data = np.array(exp_seqpos_data)


				ax.plot( pred_seqpos_data, exp_seqpos_data, marker='x', markersize=4, linestyle=' ' )
				ax.xlim(0,3)
				ax.ylim(0,3)
				ax.plot( plt.xlim(), plt.ylim(), color='black')
				plt.xlabel('Prediction (Seqpos %s)' % str(pred_seqpos), fontsize=12)
				plt.ylabel('Experimental (Seqpos %s)' % str(exp_seqpos), fontsize=12)




	# save as pdf and close
	pp.savefig()
	pp.close()

	# open pdf
	open_figure(fullpdfname)