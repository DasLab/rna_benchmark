#!/usr/bin/python

##########################################################

from os.path import exists, dirname, basename, abspath
import matplotlib.pyplot as plt
import numpy as np
from make_plots_util import *
import subprocess

##########################################################

def make_plots( inpaths, outfilenames=['swm_rebuild.out','swm_rebuild.sc'], target_files=['favorites.txt','favorites2.txt'], targets=[''], colorcode=None, xvars=['rms_fill'], yvars=['score'], scale=False ):

	data = []
	which_target = []
	outfiles_list = []

	if not colorcode: 
		colorcode = [ (0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0) ]
	if len(colorcode) < len(inpaths): 
		colorcode = jet( len(inpaths) )

	if targets[0] != '':
		target_names = targets
	else:
		target_names = get_target_names( target_files )

	for target in target_names:
		print "Target: "+target

	inpaths = map( lambda x: abspath(x), inpaths )
	base_inpaths = map( lambda x: basename(x), inpaths )

	for n in xrange( len(inpaths) ):
		assert( exists( inpaths[n] ) )
		outfiles = []
		for outfilename in outfilenames:
			outfiles += get_outfiles( inpaths[n], outfilename )
		outfiles_actual = []
		for outfile in outfiles:
			if outfile.find( '.out' ) > 0 and outfile.replace( '.out','.sc' ) in outfiles:	continue
			print 'Reading in ... '+outfile
			assert( exists( outfile ) )
			outfiles_actual.append( outfile )
	  	which_target.append( map( lambda x: target_names.index( basename( dirname( x ) ) ), outfiles_actual ) )
		data.append( dict([ (target_names[ which_target[n][k] ], load_score_data(outfiles_actual[k])) for k in xrange( len(outfiles_actual) ) ]) )
		outfiles_list.append( outfiles_actual )
	noutfiles = np.max( map( lambda x: len(x), outfiles_list ) )

	###################################################

	# get and print out runtimes, stored in the silent files
	show_times( inpaths, data, noutfiles, target_names, which_target )
	times_list = get_times( inpaths, data, noutfiles, target_names, which_target )
	
	###################################################

	# setup pdf file name, and PdfPages handle
	( pp, fullpdfname ) = setup_pdf_page( base_inpaths )
	
	# get nplots, nrows, ncols, figwidth, figheight
	( nplots, nrows, ncols, figwidth, figheight ) = get_figure_dimensions( noutfiles )

	# set up figure, adjust properties
	fig = plt.figure(1)
	fig.set_size_inches( figwidth, figheight )

	# iterate over runs
	for n in xrange( len(inpaths) ):

		# iterate over targets
		for plot_idx, target in enumerate(target_names, start=1):

			# get run_time from times_list
			run_time = times_list[n][plot_idx-1]

			# add subplot, if scores are available
			try: 
				data[n][target]
			except:	
				continue
			ax = fig.add_subplot( nrows, ncols, plot_idx )

			# get data
			for xvar in xvars:
				if xvar not in data[n][target].score_labels: continue
				xvar_idx = data[n][target].score_labels.index( xvar )
				break
			assert( xvar_idx > -1 )

			for yvar in yvars:
				if yvar not in data[n][target].score_labels: continue
				yvar_idx = data[n][target].score_labels.index( yvar )
				break
			assert( yvar_idx > -1 )
			
			[ xvar_data, yvar_data ] = [ list(d) for d in zip( *[ ( score[xvar_idx], score[yvar_idx] ) for score in data[n][ target ].scores] ) ]

			# plot data
			if run_time.times_found():
				time_label = '%5.0f +/- %4.0f' %( run_time.mean, run_time.stdev )
			else:
				time_label = ''
			plot_label = base_inpaths[n]
			ax.plot( xvar_data, yvar_data, marker='.', markersize=4, color=colorcode[n], linestyle=' ', label=plot_label )
			ax.plot( [1 for y in plt.ylim()], plt.ylim(), color='black', linestyle=':')
			ax.plot( [2 for y in plt.ylim()], plt.ylim(), color='black')

			# set axes limits
			if not scale:	
				ax.set_xlim( 0, 16 )

			# set title and axes labels		
			ax.set_title( get_title(target), fontsize='medium', weight='bold' )
			ax.set_ylabel( string.join(yvars, ', '), fontsize=6 )
			ax.set_xlabel( string.join(xvars, ', '), fontsize=6 )

			# adjust axis properties
			for tick in ( ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels() ):	
				tick.set_fontsize(6)

			# setup times sublegends
			fontsize = 6
			ax.text( 
				0.92, 
				(0.10*len(inpaths)) - (0.015*fontsize*n), 
				time_label,
        		verticalalignment='bottom', 
        		horizontalalignment='right',
        		transform=ax.transAxes,
        		color=colorcode[n], fontsize=fontsize
        	)

			# setup global legend based on inpaths
			( handles, labels ) = ax.get_legend_handles_labels()
			if (plot_idx == 1 or nplots < 3):
				legend = ax.legend(handles, base_inpaths[:n+1], loc=1, numpoints=1, prop={'size':6})

	# adjust spacing of plots on figure
	if ( nplots == 1 or nrows < ncols ): # landscape
		plt.subplots_adjust(bottom=.1, left=.05, right=.98, top=.90, hspace=.5)
	else:			
		plt.subplots_adjust(bottom=.075, left=.08, right=.95, top=.95, wspace=.3, hspace=.5)

	# print date to figure ( bottom_right = (0.99, 0.01); top_right = (0.99, 0.98) )
	plt.figtext(0.95, 0.02, get_date(), horizontalalignment='right') 

	# save as pdf and close
	pp.savefig()
	pp.close()

	# open pdf 
	try:
		subprocess.call( ['open',fullpdfname] ) # works nicely on a mac.
	except:
		pass

	return

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='Make plots of scores from silent files.')
	parser.add_argument('inpaths', nargs='+', help='List of paths too silent files.')
	parser.add_argument('-outfilenames', nargs='*', help='Name of silent file.', default=['swm_rebuild.out','swm_rebuild.sc'])
	parser.add_argument('-target_files', nargs='+', help='List of additional target files.', default=['favorites.txt','favorites2.txt'])
	parser.add_argument('-targets', nargs='+', help='List of targets.', default=[''])
	parser.add_argument('-xvar', nargs='*', help='Name of x variable(s).', default=['rms_fill'])
	parser.add_argument('-yvar', nargs='*', help='Name of y variable(s).', default=['score'])
	parser.add_argument('--scale', help='scale plot axes.', action='store_true')
	args=parser.parse_args()

	make_plots( args.inpaths, outfilenames=args.outfilenames, target_files=args.target_files, targets=args.targets, xvars=args.xvar, yvars=args.yvar, scale=args.scale )

