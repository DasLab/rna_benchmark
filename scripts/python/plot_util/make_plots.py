#!/usr/bin/python

###############################################################################
### imports
###############################################################################
import sys
import os
from os.path import exists, dirname, basename, abspath, isdir
import matplotlib.pyplot as plt
import numpy as np
from make_plots_util import *
import subprocess
from matplotlib.font_manager import FontProperties
import argparse

###############################################################################
### main functions
###############################################################################
def make_plots(argv):


	# initialize options
	options = init_options(argv)
	inpaths = options.inpaths
	outfilenames = options.outfilenames
	target_files = options.target_files
	targets = options.targets
	xvars, yvars = options.xvar, options.yvar
	pdfname = options.pdfname

        # option-dependent imports
        if options.seaborn is True:
                print "[option-dependent imports]"
                print " * options.seaborn:", options.seaborn
                print " * importing seaborn as sns"
                import seaborn as sns

	# check options
	inpaths = [abspath(x) for x in inpaths if exists(x) and isdir(x)]
	if targets[0] != '*':
		targets = targets
	elif target_files is not None:
		targets = get_target_names( target_files )
	else:
		targets = get_target_names( target_files, inpaths )
	for target in targets:
		print "Target: "+target

	# Load data for all targets in all inpaths
	data = load_data( inpaths, targets, outfilenames )

	inpaths = [x for x in inpaths if x in data.keys()]
	nplots = len(set([k for v in data.values() for k in v.keys()]))
        if ( nplots == 0 ):
                print "You need to specify -target_files, or your inpaths are nt valid."
                exit

	# get and print out runtimes, stored in the silent files
	times_list = get_times( inpaths, data, targets, verbose=True )

	###################################################

	# setup pdf and figure, return handles
	( pp, fullpdfname ) = setup_pdf_page( inpaths, targets, pdfname=pdfname )
	( fig, nplots, nrows, ncols ) = setup_figure( nplots )
	colorcode = get_colorcode( len(inpaths), options=options )
        markersize = 4

	xlabels = []
	ylabels = []
	plot_idx = 0
	handles = []
	# iterate over targets
	for target_idx, target in enumerate(targets, start=1):

		# make sure target is in data for atleast one inpath
		if not any ( target in d.keys() for d in data.values() ):
			continue

		# get subplot, if data exists for target
		plot_idx += 1
		ax = fig.add_subplot( nrows, ncols, plot_idx )

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

			# check for sequence recovery
			if "sequence_recovery" in xvars+yvars:
				score_labels.append("sequence_recovery")

				tag_idx = score_labels.index("description")
				for idx, score in enumerate(data[inpath][target].scores):

					tag = score[tag_idx]
					seq_recovery = get_sequence_recovery(inpath, target, outfilenames, tag)
					data[inpath][target].scores[idx].append(seq_recovery)

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
			for score in data[inpath][target].scores:

				if options.ignore_missing:
					missing_idx = score_labels.index("missing")
					if missing_idx > -1 and float(score[missing_idx]) > 0.0:
						continue

				xvar_data.append( float(score[xvar_idx]) )
				yvar_data.append( float(score[yvar_idx]) )

                        label = basename(inpath)
                      	if options.seaborn is True:
                                sns.set_style("darkgrid")
                                sns.set_context("poster")
                                markersize = 10

                                label = 'StepWise Assembly'# (SWA)'
                                if 'swm' in basename(inpath):
                                        label = 'StepWise Monte Carlo'# (SWM)'
                                if 'farfar' in basename(inpath).lower():
                                        label = 'FARFAR'# (SWM)'


			# plot data, and reference lines (x=1, x=2)
			h = ax.plot( xvar_data,
				 yvar_data,
				 marker='o',
				 markersize=markersize,
				 color=colorcode[inpath_idx],
				 linestyle=' ',
				 label=label )
			handles.append(h)

                        if options.seaborn is False:
                                ax.plot( [1 for y in plt.ylim()],
                                         plt.ylim(),
                                         color='black',
                                         linestyle=':')
                                ax.plot( [2 for y in plt.ylim()],
                                         plt.ylim(),
                                         color='black')

			# set title and axes labels, adjust axis properties
			title_fontsize = 'small' if nrows < 3 else 8
                        ax.set_title( get_title(target), fontsize=title_fontsize, weight='bold' )
                        ax.set_ylabel(string.join(ylabels, ', '), fontsize=6 )
                        ax.set_xlabel(string.join(xlabels, ', '), fontsize=6 )
                        for ticklabel in ax.yaxis.get_ticklabels()+ax.xaxis.get_ticklabels():
                                ticklabel.set_fontsize(6)
			ax.set_xlim(0, 16)

                        if options.seaborn is True:
                                ax.set_title( get_title(target), fontsize=20, weight='bold' )
                                ax.set_ylabel('Rosetta Energy', fontsize=20, weight='bold')
                                ax.set_xlabel(r'RMSD ($\AA$)', fontsize=20, weight='bold')
                                for ticklabel in ax.yaxis.get_ticklabels()+ax.xaxis.get_ticklabels():
                                        ticklabel.set_fontsize(14)


			# print times in plots (if available)
                        if options.seaborn is False and times_list[inpath_idx][target_idx-1].times_found():
                                monospace_font = FontProperties()
                                monospace_font.set_family( 'monospace' )
                                xpos, ypos = 0.92, (0.10*len(inpaths)) - (0.015*6*inpath_idx)
				ax.text( xpos,
					 ypos,
					 times_list[inpath_idx][target_idx-1].get_label(),
					 verticalalignment='bottom',
					 horizontalalignment='right',
					 transform=ax.transAxes,
					 color=colorcode[inpath_idx],
					 fontsize=6,
					 fontproperties=monospace_font )

	# finalize (adjust spacing, print date)
	finalize_figure( fig, nplots, nrows, ncols, options=options )

	# save as pdf and close
	pp.savefig()
	pp.close()

	# open pdf
	out, err = subprocess.Popen(['uname'], stdout=subprocess.PIPE).communicate()
	if 'Darwin' in out:
		subprocess.call(['open',fullpdfname])
	if 'Linux' in out:
		subprocess.call(['xdg-open',fullpdfname])

	return True


###############################################################################
### initialization functions
###############################################################################
def init_options(argv):
	if isinstance(argv, argparse.ArgumentParser):
		return argv
	options = init_options_parser().parse_args(args = argv)
	return options

def init_options_parser():
	parser = argparse.ArgumentParser(description='Plot scores from silent files.')
	parser.add_argument(
		'inpaths',
		nargs='+',
		help='List of paths to silent files.'
	)
	parser.add_argument(
		'-outfilenames',
		nargs='*',
		help='Name of silent file.',
		default=['swm_rebuild.out',
			 'swm_rebuild.sc',
                         'farna_rebuild.sc',
                         'farna_rebuild.out',
			 'region_FINAL.out']
	)
	parser.add_argument(
		'-target_files',
		nargs='+',
		help='List of additional target files.',
		default=['favorites.txt','favorites2.txt','challenges.txt','RNA_loop_motifs_PS2011.txt','followups.txt']
	)
	parser.add_argument(
		'-targets',
		nargs='+',
		help='List of targets.',
		default=['*']
	)
	parser.add_argument(
		'-xvar',
		nargs='*',
		help='Name of x variable(s).',
		default=['rms_fill','rms','NAT_rmsd']
	)
	parser.add_argument(
		'-yvar',
		nargs='*',
		help='Name of y variable(s).',
		default=['score']
	)
	parser.add_argument(
		'-o','--pdfname',
		help='File name to save as pdf.',
		default=None
	)
	parser.add_argument(
		'--ignore_missing',
		help='Ignore structures with unbuilt residues.',
		action = "store_true"
	)
	parser.add_argument(
		'--seaborn',
		help='Use seaborn plotting styles.',
		action = "store_true"
	)
	return parser


###############################################################################
### main
###############################################################################
if __name__=='__main__':
	sys.exit(make_plots(sys.argv[1:]))

