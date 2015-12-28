#!/usr/bin/python

##########################################################

from os.path import exists, dirname, basename, abspath
import matplotlib.pyplot as plt
import numpy as np
from make_plots_util import *
import subprocess
from matplotlib.font_manager import FontProperties

##########################################################

def make_plots( inpaths, outfilenames=[], target_files=['favorites.txt','favorites2.txt','challenges.txt'], targets=['*'], colorcode=None, xvars=['rms_fill'], yvars=['score'] ):

	# initialize lists
	data = []
	which_target = []
	outfiles_list = []
	bad_inpaths = []

	# Get absolute paths of inpaths
	inpaths = map( lambda x: abspath(x), inpaths )

	# Get colorcode for plotting
	if not colorcode:
		colorcode = [ (0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0), (0.0, 0.5, 1.0, 1.0), (0.2, 0.7, 0.0, 1.0), (1.0, 0.0, 1.0)  ]
	if len(colorcode) < len(inpaths):
		colorcode = jet( len(inpaths) )

	# Get target_names
	if targets[0] != '*':
		target_names = targets
	else:
		target_names = get_target_names( target_files )

	# Print target_names
	for target in target_names:
		print "Target: "+target

	# Check all inpaths for outfiles, get outfiles/data
	for n, inpath in enumerate(inpaths):
		assert( exists( inpath ) )
		outfiles_actual = get_outfiles( inpath, target_names, outfilenames )
		if not len(outfiles_actual):  bad_inpaths.append( inpath )
	  	which_target.append( map( lambda x: target_names.index( basename( dirname( x ) ) ), outfiles_actual ) )
		data.append( dict([ (target_names[ which_target[n][k] ], load_score_data(outfiles_actual[k])) for k in xrange( len(outfiles_actual) ) ]) )
		outfiles_list.append( outfiles_actual )

	# Remove empty items from list
	which_target = [item for item in which_target if len(item)]
	data = [item for item in data if len(item)]
	outfiles_list = [item for item in outfiles_list if len(item)]

	# Remove inpath from inpaths if no outfile found for targets in inpath
	inpaths = [ path for path in inpaths if path not in bad_inpaths ]
	base_inpaths = map( lambda x: basename(x), inpaths )

	# Get the max number of outfiles to be plotted for a single inpath
        if len( outfiles_list ) == 0:
                print 'No .out or .sc files found!'
                exit()
	noutfiles = np.max( map( lambda x: len(x), outfiles_list ) )

	###################################################

	# get and print out runtimes, stored in the silent files
	times_list = get_times( inpaths, data, noutfiles, target_names, which_target, verbose=True )

	# setup pdf and figure, return handles
	( pp, fullpdfname ) = setup_pdf_page( base_inpaths, targets )
	( fig, nplots, nrows, ncols ) = setup_figure( noutfiles )
        title_fontsize = get_title_fontsize( ncols, nplots )

	# iterate over runs
	for n in xrange( len(inpaths) ):

		# iterate over targets
		for plot_idx, target in enumerate(target_names, start=1):

			# get subplot, if data exists for target
			if target not in data[n].keys(): continue
                        plot_idx_wrapped = 1 + (plot_idx - 1) % (nrows*ncols)
			ax = fig.add_subplot( nrows, ncols, plot_idx_wrapped )

			# get index of first xvar/yvar found in score_labels
			score_labels = data[n][target].score_labels
			for xvar in xvars:
				xvar_idx = score_labels.index( xvar ) if xvar in score_labels else -1
				if xvar_idx > -1:  break
			for yvar in yvars:
				yvar_idx = score_labels.index( yvar ) if yvar in score_labels else -1
				if yvar_idx > -1:  break

			# get data from scores using xvar_idx and yvar_idx
			assert( xvar_idx > -1 and yvar_idx > -1 )
			[ xvar_data, yvar_data ] = [ list(d) for d in zip( *[ ( score[xvar_idx], score[yvar_idx] ) for score in data[n][ target ].scores] ) ]

			# plot data, and reference lines (x=1, x=2)
			ax.plot( xvar_data, yvar_data, marker='.', markersize=4, color=colorcode[n], linestyle=' ', label=base_inpaths[n] )
                        #ylim = plt.ylim()
                        #ylim = ( max( ylim[0], -500 ) , min( ylim[1],  100 ) )
                        #ax.set_ylim( ylim )
			ax.plot( [1 for y in plt.ylim()], plt.ylim(), color='black', linestyle=':')
			ax.plot( [2 for y in plt.ylim()], plt.ylim(), color='black')
			ax.set_xlim( 0, 16 )

			# set title and axes labels, adjust axis properties
			ax.set_title( get_title( target ), fontsize=title_fontsize, weight='bold' )
			ax.set_ylabel( string.join(yvars, ', '), fontsize=6 )
			ax.set_xlabel( string.join(xvars, ', '), fontsize=6 )
			for ticklabel in ax.yaxis.get_ticklabels()+ax.xaxis.get_ticklabels():
				ticklabel.set_fontsize(6)

			# print times in plots (if available)
                        monospace_font = FontProperties()
                        monospace_font.set_family( 'monospace' )
			if times_list[n][plot_idx-1].times_found():
				xpos, ypos = 0.92, (0.10*len(inpaths)) - (0.015*6*n)
				ax.text( xpos, ypos, times_list[n][plot_idx-1].get_label(),
						 verticalalignment='bottom', horizontalalignment='right',
						 transform=ax.transAxes, color=colorcode[n], fontsize=6,
                                                 fontproperties = monospace_font )

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

	parser = argparse.ArgumentParser(description='Make plots of scores from silent files.')
	parser.add_argument('inpaths', nargs='+', help='List of paths too silent files.')
	parser.add_argument('-outfilenames', nargs='*', help='Name of silent file.', default=['swm_rebuild.out','swm_rebuild.sc','farna_rebuild.sc','farna_rebuild.out'])
	parser.add_argument('-target_files', nargs='+', help='List of additional target files.', default=['favorites.txt','favorites2.txt','challenges.txt'])
	parser.add_argument('-targets', nargs='+', help='List of targets.', default=['*'])
	parser.add_argument('-xvar', nargs='*', help='Name of x variable(s).', default=['rms_fill','rms'])
	parser.add_argument('-yvar', nargs='*', help='Name of y variable(s).', default=['score'])
	args=parser.parse_args()

	make_plots( args.inpaths, outfilenames=args.outfilenames, target_files=args.target_files, targets=args.targets, xvars=args.xvar, yvars=args.yvar )
