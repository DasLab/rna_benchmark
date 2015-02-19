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
parser.add_argument('-seqfile', help='.seq file', default=None)
parser.add_argument('-exp_rdat_files', nargs='*', help='RDAT files for processing', default=[])
parser.add_argument('-pred_rdat_files', nargs='*', help='RDAT files for processing', default=[])
parser.add_argument('-exp_seqpos', nargs='*', help='list of sequence positions to compare', default=[])
parser.add_argument('-pred_seqpos', nargs='*', help='list of sequence positions to compare', default=[])

args = parser.parse_args()

###############################################################################

rdat = RDATFile()


seqfile = open( args.seqfile ).readlines()

#experimental_data = {}
#prediction_data = {}

#for rdat_file in args.exp_rdat_files:
#	experimental_data[ rdat_file ] = rdat.get_reactivity(rdat_file)

#for rdat_file in args.pred_rdat_files:
#	prediction_data[ rdat_file ] = rdat.get_reactivity(rdat_file)

#exp_seqpos_list = map(lambda x: int(x), args.exp_seqpos)
#pred_seqpos_list = map(lambda x: int(x), args.pred_seqpos)

pdfname = 'flanking_bp_heat_maps.pdf'
title = basename( pdfname.replace('_',' ').replace('.pdf','') )

########################################################################################
#flanking_bp_heat_maps( exp_seqfile, pred_seqfile, experimental_data, prediction_data, exp_seqpos, pred_seqpos, fullpdfname=pdfname, title=title )
########################################################################################
########################################################################################
#def flanking_bp_heat_maps( exp_seqfile, pred_seqfile, experimental_data, prediction_data, exp_seqpos_list, pred_seqpos_list, fullpdfname=None, title=None ):
########################################################################################
fullpdfname = pdfname
print '\nMaking figure in: %s\n' % fullpdfname
pp = PdfPages( fullpdfname )

fig = plt.figure()
fig.set_size_inches(8.5, 11)

nrows = 2
ncols = 3

# get flanking base pairs of each sample
flanking_bps = {}
for idx, line in enumerate(seqfile,start=1):
	line = line.split()
	name = line[0].split('-')
	seq1 = name[-2]
	seq2 = name[-1]
	fbp1 = '%s-%s' %(seq1[0],seq2[-1])
	fbp2 = '%s-%s' %(seq1[-1],seq2[0])
	flanking_bps[ idx ] = ( fbp1, fbp2 )

plot_idx = 0
for file_idx, rdat_file in enumerate( (args.exp_rdat_files + args.pred_rdat_files), start=1):

	print 'rdat_file = ', rdat_file
	reactivity_data = rdat.get_reactivity(rdat_file)

	if len(rdat.sequence) > 20:
		seqpos_list = [ 36, 50 ]
	else:
		seqpos_list = [ 5, 13 ]

	for idx, seqpos in enumerate(seqpos_list,start=1):

		plot_idx += 1
		#print plot_idx
		# VARIABLES TO BE DETERMINE
		#for sample_idx in flanking_bps.keys():
		#	print 'SAMPLE_IDX: ', sample_idx 
		#	print 'FLANKING_BPS_1: ', flanking_bps[ sample_idx ][0]
		#	print 'FLANKING_BPS_2: ', flanking_bps[ sample_idx ][1]

		x_bps_axis = [ 'u-a', 'a-u', 'g-c', 'c-g', 'u-g', 'g-u' ]
		y_bps_axis = [ 'u-g', 'g-u', 'g-c', 'c-g', 'u-a', 'a-u' ] #list(reversed(x_bps_axis))
		data = [[0.0 for x in x_bps_axis] for y in y_bps_axis]
		#print data

		for sample_idx in flanking_bps.keys():
			yidx = y_bps_axis.index( flanking_bps[ sample_idx ][0] )
			xidx = x_bps_axis.index( flanking_bps[ sample_idx ][1] )
			data[ yidx ][ xidx ] += reactivity_data[ sample_idx ][ seqpos ]

		data = np.array( data )
		#print data 

		ax = fig.add_subplot( ncols, nrows, plot_idx )
		plt.pcolor(data, cmap=plt.get_cmap('Greys'))
		plt.xticks(np.arange(0,len(x_bps_axis))+0.5)
		plt.yticks(np.arange(0,len(y_bps_axis))+0.5)
		#ax.xaxis.tick_top()
		#ax.yaxis.tick_left()

		ax.set_xticklabels( x_bps_axis, minor=False, fontsize=12)
		ax.set_yticklabels( y_bps_axis, minor=False, fontsize=12)
		ax.set_aspect( 1 )

		ax.set_xlabel('Bottom Base Pair')
		ax.set_ylabel('Top Base Pair')

		title = rdat_file.replace('.rdat','').replace('_',' ')
		if idx == 1:
			title += " (Seqpos 5'-A)"
		else:
			title += " (Seqpos 3'-A)"

		ax.set_title(title, fontsize=8, weight='bold')

		plt.subplots_adjust(hspace=0.35, wspace=0.35)



# save as pdf and close
pp.savefig()
pp.close()

# open pdf
open_figure(fullpdfname)