#!/usr/bin/python

##########################################################

from sys import exit
from os.path import exists, dirname, basename, abspath
from os import popen
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from make_plots_util import *

##########################################################

def make_plots( inpaths, outfilename='swm_rebuild.out', target_files=['favorites.txt','favorites2.txt'], colorcode=None, xvar='rms_fill', yvar='score', scale=False, show=False, landscape=False ):
	
	data = []
	which_target = []
	outfiles_list = []
	
	if not colorcode: colorcode = [ (0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0) ]
	if len( colorcode ) < len( inpaths ): colorcode = jet( len( inpaths ) )
	target_names = get_target_names( target_files )
	inpaths = map( lambda x: abspath(x), inpaths )
	
	for n in xrange( len(inpaths) ):
		assert( exists( inpaths[n] ) )
		outfiles = popen( 'ls -1 '+inpaths[n]+'/*/'+outfilename ).read().split('\n')[:-1]
		for outfile in outfiles:
			print 'Reading in ... '+outfile
			assert( exists( outfile ) )
		which_target.append( map( lambda x: target_names.index( basename( dirname( x ) ) ), outfiles ) )
		data.append( dict([ (target_names[ which_target[n][k] ], load_score_data(outfiles[k])) for k in xrange( len(outfiles) ) ]) )
		outfiles_list.append( outfiles )
	noutfiles = np.max( map( lambda x: len(x), outfiles_list ) )

	###################################################
	
	# print out runtimes, stored in the silent files
	show_times( inpaths, data, noutfiles, target_names, which_target )

	###################################################

	# get nplots, nrows, ncols
	nplots = noutfiles
	assert( nplots )
	if landscape:
		if nplots < 3: 	  nrows = 1
		else: 			  nrows = 3
	else:
		if nplots < 3: 	  nrows = nplots
		elif nplots < 12: nrows = 4
		else:  		      nrows = 5
	ncols = np.ceil( nplots / float( nrows ) )

	# setup pdf file name, and PdfPages handle
	pdfname = basename( inpaths[0] )
	if len( inpaths ) > 1:	
		for k in xrange( 1, len( inpaths ) ): pdfname += '_vs_' + basename( inpaths[k] )
	fullpdfname = get_path_to_dir('stepwise_benchmark') + '/Figures/' + pdfname # + '.pdf'
	if landscape:	fullpdfname += '_landscape.pdf'
	else:			fullpdfname += '.pdf'
	print '\nMaking figure in: %s\n' % fullpdfname
	pp = PdfPages( fullpdfname )

	# set up figure, adjust properties
	fig = plt.figure(1)
	if landscape:	fig.set_size_inches(11,8.5)
	else:			fig.set_size_inches(8.5,11)

	# iterate over runs
	for n in xrange( len( inpaths ) ):

		# initialize plot index
		plot_idx = 0

		# iterate over targets
		for target in target_names:
			
			# add subplot, if scores are available
			plot_idx += 1
			try: data[n][ target ]
			except:	continue
			ax = fig.add_subplot( nrows, ncols, plot_idx )

			# get data
			( xvar_idx , yvar_idx  ) = data[n][ target ].score_labels.index( xvar ) , data[n][ target ].score_labels.index( yvar )
			[ xvar_data, yvar_data ] = [ list(d) for d in zip( *[ ( score[xvar_idx], score[yvar_idx] ) for score in data[n][ target ].scores] ) ]

			# plot data
			ax.plot( xvar_data, yvar_data, marker='.', markersize=5, color=colorcode[n], linestyle=' ', label=basename(inpaths[n]) )	
			ax.plot( [1 for y in plt.ylim()], plt.ylim(), color='black', linestyle=':')
			ax.plot( [2 for y in plt.ylim()], plt.ylim(), color='black')

			# set axes limits
			if not scale:	ax.set_xlim( 0, 16 )

			# set title and axes lables
			if landscape:	ax.set_title( get_title(target), fontsize='small', weight='bold' )
			else:			ax.set_title( get_title(target), fontsize='medium', weight='bold' )
			if ( ( np.mod( plot_idx, ncols ) == 1 ) or ( ncols == 1 ) ):	
				if landscape:	ax.set_ylabel( yvar, fontsize='small' )
				else:			ax.set_ylabel( yvar, fontsize='medium' )
			if ( ( np.floor( (plot_idx-1) / ncols ) == nrows-1 ) or ( nrows == 1 ) ):	
				if landscape:	ax.set_xlabel( xvar, fontsize='small' )
				else:			ax.set_xlabel( xvar, fontsize='medium' )

			# adjust axis properties
			for tick in ax.xaxis.get_ticklabels():	tick.set_fontsize(6)
			for tick in ax.yaxis.get_ticklabels():	tick.set_fontsize(6)

			# setup legend
			if ( plot_idx == 3 ):	
				legend = ax.legend(shadow=True)
				for label in legend.get_texts():	label.set_fontsize(8)
				for label in legend.get_lines():	label.set_linewidth(.5)

	# adjust spacing of plots on figure
	if landscape:	plt.subplots_adjust(bottom=.1, left=.05, right=.98, top=.90, hspace=.5)
	else:			plt.subplots_adjust(bottom=.05, left=.08, right=.95, top=.95, hspace=.35)

	# save as pdf and close, show plot if show=True
	pp.savefig()
	pp.close()
	if show:	plt.show()

	return

##########################################################
##########################################################

if __name__=='__main__':

	import argparse

	parser = argparse.ArgumentParser(description='Make plots of scores from silent files.')
	parser.add_argument('inpaths', nargs='+', help='List of paths to silent files.')
	parser.add_argument('-outfilename', help='Name of silent file.', default='swm_rebuild.out')
	parser.add_argument('-target_files', nargs='+', help='List of additional target files.', default=['favorites.txt','favorites2.txt'])
	parser.add_argument('-xvar', help='Name of x variable.', default='rms_fill')
	parser.add_argument('-yvar', help='Name of y variable.', default='score')
	parser.add_argument('--scale', help='scale plot axes.', action='store_true')
	parser.add_argument('--landscape', help='orientation of figure.', action='store_true')
	parser.add_argument('--show', help='show plot after it is created.', action='store_true')
	args=parser.parse_args()

	make_plots( args.inpaths, outfilename=args.outfilename, target_files=args.target_files, xvar=args.xvar, yvar=args.yvar, scale=args.scale, show=args.show, landscape=args.landscape )
