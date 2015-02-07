#!/usr/local/bin/python

###############################################################################

import os
import subprocess
from os.path import basename, dirname, exists, isfile, isdir
import numpy as np

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
				continue
			if 'NAME' in line:
				self.name = line.split()[1:]
				continue
			if 'SEQUENCE' in line:
				self.sequence = line.split()[1:]
				continue
			if 'STRUCTURE' in line:
				self.structure = line.split()[1:]
				continue
			if 'ANNOTATION' in line:
				#self.annotation = line.split()[1:]
				continue
			if 'SEQPOS' in line:
				seqpos_list = line.split()[1:]
				for idx, seqpos in enumerate(seqpos_list, start=1):
					self.seqpos[ idx ] = seqpos.lower()
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


	def total_seqpos_reactivity(self, rdat_files, slice_seqpos=None):
		total_seqpos_reactivity = {}
		for rdat_file in rdat_files:
			self.load(rdat_file)
			for sample_idx in sorted(self.reactivity.keys()):
				if not sample_idx in total_seqpos_reactivity.keys():
					total_seqpos_reactivity[ sample_idx ] = {}
				for seqpos_idx in sorted(self.reactivity[ sample_idx ].keys()):
					if slice_seqpos:
						if seqpos_idx not in slice_seqpos:
							continue 
					seqpos_reactivity = self.reactivity[ sample_idx ][ seqpos_idx ]
					if not seqpos_idx in total_seqpos_reactivity[ sample_idx ].keys():
						total_seqpos_reactivity[ sample_idx ][ seqpos_idx ] = seqpos_reactivity
					else:
						total_seqpos_reactivity[ sample_idx ][ seqpos_idx ] += seqpos_reactivity
		return total_seqpos_reactivity

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
		else:
			return mean_seqpos_reactivity

	def seqpos_reactivity(self, rdat_file, slice_seqpos=None, verbose=False):
		seqpos_reactivity = self.total_seqpos_reactivity([rdat_file], slice_seqpos=slice_seqpos)				
		for sample_idx in sorted(seqpos_reactivity.keys()):
			if verbose:
				print 'SAMPLE:', sample_idx
				print 'SEQPOS_IDX  MEAN_REACTIVITY'
				for seqpos_idx in list(reversed(sorted(seqpos_reactivity[ sample_idx ].keys()))):
					print '%-10d  %f' %( seqpos_idx, seqpos_reactivity[ sample_idx ][ seqpos_idx ] )
		return seqpos_reactivity

	def make_rna_chemical_map(self, seqpos_reactivity_data, ndecoys=None):

		sample_idx_list = sorted(seqpos_reactivity_data.keys())
		seqpos_idx_list = sorted(seqpos_reactivity_data[ sample_idx_list[0] ].keys())	
		
		data = []
		for seqpos_idx in seqpos_idx_list:
			seqpos_data = [ seqpos_reactivity_data[ x ][ seqpos_idx ] for x in sample_idx_list ]
			data.append(np.array(seqpos_data)) 
		data = np.array(data)

		rows = seqpos_idx_list
		columns = [x for x in xrange(1, len(seqpos_reactivity_data.keys())+1)]

		print 'DATA: ', data
		print 'ROWS: ', rows
		print 'COLS: ', columns

		save_fig_dir = '../../../Figures/'
		run_dir = basename(os.getcwd())
		fullpdfname = '%s/DMS_reactivity_predictions_%s.pdf' % (save_fig_dir, run_dir)
		if ndecoys:
			fullpdfname.replace('.pdf', '_%d_decoys.pdf'%ndecoys)
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


if __name__=='__main__':

	import sys

	try:
		rdat_file = sys.argv[1]
	except:
		print 'invalid argv[1]'
		sys.exit(-1)

	if len(sys.argv) > 2:
		start1 = int(sys.argv[2])
		stop1 = int(sys.argv[3])
		slice_seqpos = [x for x in xrange(start1,stop1+1)]
		if len(sys.argv) > 5:
			start2 = int(sys.argv[4])
			stop2 = int(sys.argv[5])
			slice_seqpos += [x for x in xrange(start2,stop2+1)]

	else:
		seqpos_range = None
	rdat = RDATFile()

	target_seqpos_reactivity = rdat.seqpos_reactivity(rdat_file, slice_seqpos=slice_seqpos, verbose=True)
	rdat.make_rna_chemical_map(target_seqpos_reactivity)




