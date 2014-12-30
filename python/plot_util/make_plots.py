#!/usr/bin/python

##########################################################

from os.path import exists, dirname, basename, abspath
import matplotlib.pyplot as plt
import numpy as np
from make_plots_util import *
import subprocess

##########################################################

def make_plots( inpaths, outfilenames=['swm_rebuild.out'], target_files=['favorites.txt','favorites2.txt'], targets=[''], colorcode=None, xvar='rms_fill', yvar='score', scale=False, show=False, landscape=False ):

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
		outfiles = []
		outfiles_actual = []
		assert( exists( inpaths[n] ) )
    	for outfilename in outfilenames:
        	outfiles += get_outfiles( inpaths[n], outfilename )
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
	pp = setup_pdf_page( base_inpaths, landscape=landscape )
	
	# get nplots, nrows, ncols, figwidth, figheight
	( nplots, nrows, ncols, figwidth, figheight ) = get_figure_dimensions( noutfiles, landscape=landscape )

	# set up figure, adjust properties
	fig = plt.figure(1)
	fig.set_size_inches( figwidth, figheight )

	# iterate over runs
	for n in xrange( len(inpaths) ):

		# initialize plot index
		plot_idx = 0

		# iterate over targets
		for target in target_names:

			# get run_time from times_list
			run_time = times_list[n][plot_idx]

			# add subplot, if scores are available
			plot_idx += 1
			try: 
				data[n][target]
			except:	
				continue
			ax = fig.add_subplot( nrows, ncols, plot_idx )

			# get data
			( xvar_idx , yvar_idx  ) = data[n][ target ].score_labels.index( xvar ) , data[n][ target ].score_labels.index( yvar )
			[ xvar_data, yvar_data ] = [ list(d) for d in zip( *[ ( score[xvar_idx], score[yvar_idx] ) for score in data[n][ target ].scores] ) ]

			# plot data
			time_label = '%5.0f +/- %4.0f' %( run_time.mean, run_time.stdev ) 
			plot_label = time_label
			ax.plot( xvar_data, yvar_data, marker='.', markersize=4, color=colorcode[n], linestyle=' ', label=plot_label )
			ax.plot( [1 for y in plt.ylim()], plt.ylim(), color='black', linestyle=':')
			ax.plot( [2 for y in plt.ylim()], plt.ylim(), color='black')

			# set axes limits
			if not scale:	
				ax.set_xlim( 0, 16 )

			# set title and axes labels
			if landscape:	
				ax.set_title( get_title(target), fontsize='small', weight='bold' )
			else:			
				ax.set_title( get_title(target), fontsize='medium', weight='bold' )
			ax.set_ylabel( yvar, fontsize=6 )
			ax.set_xlabel( xvar, fontsize=6 )

			# adjust axis properties
			for tick in ax.xaxis.get_ticklabels():	
				tick.set_fontsize(6)
			for tick in ax.yaxis.get_ticklabels():	
				tick.set_fontsize(6)

			# setup times sublegends
			#time_legend = ax.legend(loc=4, numpoints=1, prop={'size':6} )
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
				#plt.gca().add_artist(time_legend)

	# adjust spacing of plots on figure
	if landscape:	
		plt.subplots_adjust(bottom=.1, left=.05, right=.98, top=.90, hspace=.5)
	else:			
		plt.subplots_adjust(bottom=.05, left=.08, right=.95, top=.95, wspace=.3, hspace=.6)

	# print date to figure ( bottom_right = (0.99, 0.01); top_right = (0.99, 0.98) )
	plt.figtext(0.99, 0.01, get_date(), horizontalalignment='right') 

	# save as pdf and close, show plot if show=True
	pp.savefig()
	pp.close()
	if show:	
		plt.show()

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
	parser.add_argument('-xvar', help='Name of x variable.', default='rms_fill')
	parser.add_argument('-yvar', help='Name of y variable.', default='score')
	parser.add_argument('--scale', help='scale plot axes.', action='store_true')
	parser.add_argument('--landscape', help='orientation of figure.', action='store_true')
	parser.add_argument('--show', help='show plot after it is created.', action='store_true')
	args=parser.parse_args()

	make_plots( args.inpaths, outfilenames=args.outfilenames, target_files=args.target_files, targets=args.targets, xvar=args.xvar, yvar=args.yvar, scale=args.scale, show=args.show, landscape=args.landscape )

