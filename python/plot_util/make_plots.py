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

def make_plots( inpaths, outfilename='swm_rebuild.out', target_files=['favorites.txt','favorites2.txt'], colorcode=None, xvar='rms_fill', yvar='score', scale=False, show=False ):

	if not colorcode: colorcode = [ (0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0) ]
	if len( colorcode ) < len( inpaths ): colorcode = jet( len( inpaths ) )

	target_names = get_target_names( target_files )

	data = []
	which_target = []
	outfiles_list = []

	inpaths = map( lambda x: abspath(x), inpaths )

	for inpath in inpaths:

		assert( exists( inpath ) )
		outfiles = popen( 'ls -1 '+inpath+'/*/'+outfilename ).read().split('\n')[:-1]

		for outfile in outfiles:
			print 'Reading in ... '+outfile
			assert( exists( outfile ) )

		data.append( map( lambda x: load_score_data( x ), outfiles ) )
		which_target.append( map( lambda x: target_names.index( basename( dirname( x ) ) ), outfiles ) )
		outfiles_list.append( outfiles )

	noutfiles = np.max( map( lambda x: len(x), outfiles_list ) )

	###################################################

	time_name = 'time'
	times_list = []

	for n in xrange( len( inpaths ) ):
		times = []
		for k in xrange( noutfiles ):
			try:
				time_idx = data[n][k].score_labels.index( time_name )
				times.append( np.array( [ score[time_idx] for score in data[n][k].scores ], dtype = 'float_' ) )
			except:
				times.append( np.array( [ 0.0 for score in xrange( noutfiles ) ], dtype = 'float_' ) )
		times_list.append( times )


	for k in xrange( noutfiles ):
		print '\n %-6s%33s' % ( 'TARGET', target_names[ which_target[0][k] ] )
		for n in xrange( len( inpaths ) ):
			mean_time = np.mean( times_list[n][k] )
			std_time = np.std( times_list[n][k] )
			print ' Run %d                    %5.0f +/- %4.0f' % ( n, mean_time, std_time )
	print '\n'
	for n in xrange( len( inpaths ) ):
		print ' Run %d: %s' % ( n, basename( inpaths[n] ) )

	###################################################

	pdfname = basename( inpaths[0] )
	if len( inpaths ) > 1:
		for k in xrange( 1, len( inpaths ) ): pdfname += '_vs_' + basename( inpaths[k] )
	fullpdfname = get_path_to_dir('stepwise_benchmark') + '/Figures/' + pdfname + '.pdf'
	print '\nMaking figure in: %s\n' % fullpdfname
	pp = PdfPages( fullpdfname )

	###################################################

	
	### determine number of plots, rows, and columns
	#nplots = 0
	#for n in xrange( len( inpaths ) ):
	#	for k in xrange( len( outfiles_list[n] ) ):
	#		if not len( data[n][k].scores ): continue
	#		nplots += 1

	nplots = noutfiles
	assert( nplots )
	if nplots < 3: 	  nrows = nplots
	elif nplots < 12: nrows = 4
	else:  		      nrows = 5
	ncols = np.ceil( nplots / float( nrows ) )

	
	fig = plt.figure(1) #figsize=(1.5*ncols, 1.5*nrows))
	fig.set_size_inches(8.5,11)

	titles = []

	for n in xrange( len( inpaths ) ):

		plot_idx = 0

		for k in xrange( len( outfiles_list[n] ) ):

			if not len( data[n][k].scores ): continue
			plot_idx += 1
			#plt.subplot( nrows, ncols, plot_idx )
			ax = fig.add_subplot( nrows, ncols, plot_idx )

			( xvar_idx , yvar_idx  ) = data[n][k].score_labels.index( xvar ) , data[n][k].score_labels.index( yvar )
			[ xvar_data, yvar_data ] = [ list(d) for d in zip( *[ ( score[xvar_idx], score[yvar_idx] ) for score in data[n][k].scores] ) ]

			ax.plot( xvar_data, yvar_data, marker='.', markersize=5, color=colorcode[n], linestyle=' ', label=basename(inpaths[n]) )
			ax.set_title( target_names[ which_target[n][k] ], fontsize='medium', weight='bold' )

			if not scale:	plt.xlim( 0, 12 )

			if ( ( np.mod( plot_idx, ncols ) == 1 ) or ( ncols == 1 ) ):
				ax.set_ylabel( yvar, fontsize='medium' )
			if ( ( np.floor( (plot_idx-1) / ncols ) == nrows-1 ) or ( nrows == 1 ) ):
				ax.set_xlabel( xvar, fontsize='medium' )

			for tick in ax.xaxis.get_ticklabels() :
				tick.set_fontsize(6)

			for tick in ax.yaxis.get_ticklabels() :
				tick.set_fontsize(6)

			if ( plot_idx == 1 ):	
				ax.legend()
				plt.rc('legend', fontsize=6)


		titles.append( basename( inpaths[n] ) )

	#plt.tight_layout()
	plt.subplots_adjust(bottom=.05, left=.08, right=.95, top=.95, hspace=.35)
	ax = fig.add_subplot( nrows, ncols, 1 )
	#ax.legend( titles, prop={'size':6} )

	pp.savefig()
	pp.close()
	if show:
		plt.show()

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
	parser.add_argument('--show', help='show plot after it is created.', action='store_true')
	args=parser.parse_args()

	make_plots( args.inpaths, outfilename=args.outfilename, target_files=args.target_files, xvar=args.xvar, yvar=args.yvar, scale=args.scale, show=args.show )
