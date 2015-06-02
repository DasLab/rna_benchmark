#!/usr/bin/python

##########################################################

from os.path import exists, dirname, basename, abspath, isdir
import matplotlib.pyplot as plt
import numpy as np
from make_plots_util import *
import subprocess
from matplotlib.font_manager import FontProperties

##########################################################

def make_plots( inpaths, outfilenames, target_files, targets, xvars, yvars, pdfname ):

	# Initialize/Check args
	inpaths = [abspath(x) for x in inpaths if exists(x) and isdir(x)]
	targets = targets if targets[0] != '*' else get_target_names( target_files )
	for target in targets:
		print "Target: "+target

	# Load data for all targets in all inpaths
	data = load_data( inpaths, targets, outfilenames )
	inpaths = [x for x in inpaths if x in data.keys()]
	nplots = len(set([k for v in data.values() for k in v.keys()]))

	# get and print out runtimes, stored in the silent files
	times_list = get_times( inpaths, data, targets, verbose=True )

	###################################################

	# setup pdf and figure, return handles
	( pp, fullpdfname ) = setup_pdf_page( inpaths, targets, pdfname=pdfname )
	( fig, nplots, nrows, ncols ) = setup_figure( nplots )
	colorcode = get_colorcode( len(inpaths) )

	xlabels = []
	ylabels = []

	# iterate over targets
	for plot_idx, target in enumerate(targets, start=1):

		# make sure target is in data for atleast one inpath
		if not any ( target in d.keys() for d in data.values() ): continue

		# get subplot, if data exists for target
		plot_idx_wrapped = 1 + (plot_idx - 1) % (nrows*ncols)
		ax = fig.add_subplot( nrows, ncols, plot_idx_wrapped )

		# iterate over runs
		for inpath_idx, inpath in enumerate(inpaths):

			# plot run, if data exists
			if target not in data[inpath].keys():
				ax.plot( None,
					 None,
					 marker='.',
					 markersize=4,
					 color=colorcode[inpath_idx],
					 linestyle=' ',
					 label=basename(inpath) )
				continue

			# get index of first xvar/yvar found in score_labels
			score_labels = data[inpath][target].score_labels
			xvar_idx = -1
			for xvar in xvars:
				if not xvar in score_labels:
					continue
				xvar_idx = score_labels.index( xvar )
				if not xvar in xlabels:
					xlabels.append( xvar )
				break
			yvar_idx = -1
			for yvar in yvars:
				if not yvar in score_labels:
					continue
				yvar_idx = score_labels.index( yvar )
				if not yvar in ylabels:
					ylabels.append( yvar )
				break

			# get data from scores using xvar_idx and yvar_idx
			assert( xvar_idx > -1 and yvar_idx > -1 )
			xvar_data, yvar_data = [], []
			#xvar_cutoff, yvar_cutoff = 100.0, 5.0
			for score in data[inpath][target].scores:
				#if ( float(score[xvar_idx]) > xvar_cutoff or
				#     float(score[yvar_idx]) > yvar_cutoff ):
				#	continue
				xvar_data.append( score[xvar_idx] )
				yvar_data.append( score[yvar_idx] )

			# plot data, and reference lines (x=1, x=2)
			ax.plot( xvar_data,
				 yvar_data,
				 marker='.',
				 markersize=3,
				 color=colorcode[inpath_idx],
				 linestyle=' ',
				 label=basename(inpath) )

			ax.plot( [1 for y in plt.ylim()],
				 plt.ylim(),
				 color='black',
				 linestyle=':')
			ax.plot( [2 for y in plt.ylim()],
				 plt.ylim(),
				 color='black')

			ax.set_xlim( 0, 16 )

			# set title and axes labels, adjust axis properties
			title_fontsize = 'small' if nrows < 3 else 8
			ax.set_title( get_title(target), fontsize=title_fontsize, weight='bold' )
			ax.set_ylabel( string.join(ylabels, ', '), fontsize=6 )
			ax.set_xlabel( string.join(xlabels, ', '), fontsize=6 )
			for ticklabel in ax.yaxis.get_ticklabels()+ax.xaxis.get_ticklabels():
				ticklabel.set_fontsize(6)

			# print times in plots (if available)
			monospace_font = FontProperties()
			monospace_font.set_family( 'monospace' )
			if False:#if times_list[inpath_idx][plot_idx-1].times_found():
				xpos, ypos = 0.92, (0.10*len(inpaths)) - (0.015*6*inpath_idx)
				ax.text( xpos,
					 ypos,
					 times_list[inpath_idx][plot_idx-1].get_label(),
					 verticalalignment='bottom',
					 horizontalalignment='right',
					 transform=ax.transAxes,
					 color=colorcode[inpath_idx],
					 fontsize=6,
                                         fontproperties=monospace_font )

	# finalize (adjust spacing, print date)
	finalize_figure( fig, nplots, nrows, ncols )

	# save as pdf and close
	pp.savefig()
	pp.close()

	# open pdf
	out, err = subprocess.Popen(['uname'], stdout=subprocess.PIPE).communicate()
	if 'Darwin' in out:
		subprocess.call(['open',fullpdfname])
	if 'Linux' in out:
		subprocess.call(['xdg-open',fullpdfname])

	return

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='Plot scores from silent files.')
	parser.add_argument('inpaths',
			    nargs='+',
			    help='List of paths too silent files.'
			    )
	parser.add_argument('-outfilenames',
			    nargs='*',
			    help='Name of silent file.',
			    default=['swm_rebuild.out',
				     'swm_rebuild.sc',
				     'region_FINAL.out']
			    )
	parser.add_argument('-target_files',
			    nargs='+',
			    help='List of additional target files.',
			    default=['favorites.txt','favorites2.txt']
			    )
	parser.add_argument('-targets',
			    nargs='+',
			    help='List of targets.',
			    default=['*'])
	parser.add_argument('-xvar',
			    nargs='*',
			    help='Name of x variable(s).',
			    default=['rms_fill','NAT_rmsd']
			    )
	parser.add_argument('-yvar',
			    nargs='*',
			    help='Name of y variable(s).',
			    default=['score']
			    )
	parser.add_argument('-o','--pdfname',
			    help='File name to save as pdf.',
			    default=None
			    )
	args=parser.parse_args()


	make_plots( args.inpaths,
		    args.outfilenames,
		    args.target_files,
		    args.targets,
		    args.xvar,
		    args.yvar,
		    args.pdfname )
