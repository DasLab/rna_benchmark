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
rdat_files = args.exp_rdat_files + args.pred_rdat_files

pdfname = 'flanking_bp_heat_maps.pdf'
title = basename( pdfname.replace('_',' ').replace('.pdf','') )

fullpdfname = pdfname
print '\nMaking figure in: %s\n' % fullpdfname
pp = PdfPages( fullpdfname )

fig = plt.figure()
fig.set_size_inches(8.5, 11)

nrows = len(args.pred_seqpos)
ncols = len(rdat_files)

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
experimental = False
for file_idx, rdat_file in enumerate( rdat_files, start=1):

	print 'rdat_file:', rdat_file
	reactivity_data = rdat.get_reactivity(rdat_file)

	# need a better way of identifying whether experimental or prediction
	# so we know which seqpos to use
	if 'MgCl2' in rdat_file or 'NaCl' in rdat_file:
		seqpos_list = args.exp_seqpos #[ 36, 50 ]
		experimental = True 
	else:
		seqpos_list = args.pred_seqpos #[ 5, 13 ]
		experimental = False 

	seqpos_list = [int(x) for x in seqpos_list]

	top_bp = (seqpos_list[0]-2, seqpos_list[1]+1)
	bot_bp = (seqpos_list[0]+1, seqpos_list[1]-2)

	rxdatalist = [ reactivity_data[x][y] for y in seqpos_list for x in reactivity_data.keys()]
	cmin = min(rxdatalist)
	cmax = max(rxdatalist)

	for idx, seqpos in enumerate(seqpos_list,start=1):

		plot_idx += 1
		print 'seqpos:', seqpos
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
			data[ yidx ][ xidx ] = reactivity_data[ sample_idx ][ seqpos ]

		data = np.array( data )
		print data 

		ax = fig.add_subplot( ncols, nrows, plot_idx )
		plt.pcolor(data, cmap=plt.get_cmap('Greys'))
		print 'MIN: ', cmin
		print 'MAX: ', cmax 
		plt.clim(cmin,cmax)
		plt.colorbar()
		plt.xticks(np.arange(0,len(x_bps_axis))+0.5)
		plt.yticks(np.arange(0,len(y_bps_axis))+0.5)
		#ax.xaxis.tick_top()
		#ax.yaxis.tick_left()

		ax.set_xticklabels( [x.upper() for x in x_bps_axis], minor=False, fontsize=8, weight='bold')
		ax.set_yticklabels( [x.upper() for x in y_bps_axis], minor=False, fontsize=8, weight='bold')
		ax.set_aspect( 1 )

		ax.set_xlabel('Base Pair (%d,%d)' % (bot_bp[0],bot_bp[1]), fontsize=10 )
		ax.set_ylabel('Base Pair (%d,%d)' % (top_bp[0],top_bp[1]), fontsize=10 )

		title = rdat_file.replace('.rdat','').replace('_',' ').replace('DMS','')
		if idx == 1:
			title += " [5'-A]"#(%d)" % seqpos
		else:
			title += " [3'-A]"#%d" % seqpos

		ax.set_title(title, fontsize=10, weight='bold')

		plt.subplots_adjust(hspace=0.35, wspace=0.35)



# save as pdf and close
pp.savefig()
pp.close()

# open pdf
open_figure(fullpdfname)