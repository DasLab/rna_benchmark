#!/usr/bin/python

##########################################################

from sys import exit
from os import system, popen
from os.path import exists, dirname, basename, abspath
import string
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np

##########################################################

class ScoreData(object):

	def __init__(self):
		self.scores = None
		self.score_labels = None

##########################################################

class TimeData(object):
	
	def __init__(self):
		self.inpath = None
		self.target = None
		self.times = []
		self.mean = None
		self.stdev = None

	def update_stats(self):
		if self.times_found():
			self.mean = np.mean( np.array( self.times ) )
			self.stdev = np.std( np.array( self.times ) )
		else:
			self.mean = None
			self.stdev = None 

	def times_found(self):
		return ( len(self.times) > 0.0 ) 

##########################################################

def get_outfiles( inpath, outfilename ):
	outfiles = popen( 'ls -1 '+inpath+'/*/'+outfilename ).read()
	if 'ls:' in outfiles: # catch bad calls
		outfiles = [] 
	else: 
		outfiles = outfiles.split('\n')[:-1]
	return outfiles

##########################################################

def load_score_data( file ):
	data = ScoreData()
	for line in open( file, 'r' ).readlines():
		if len( line ) < 2: return []
		cols = line.split()
		if len( cols ) > 0 and cols[0] == 'SCORE:':
			score_labels = cols[1:len( cols )-1]
			break
	keyword = 'SCORE'
	scorenums = [ str(x) for x in xrange( 2, len( score_labels )+2 ) ] 
	command = (	
			'grep '+keyword+' '+file
			+ ' | awk \'{print $'+string.join( scorenums, ',$' )+'}\''
			+ ' | grep -v inp'	
			+ ' | grep -v R'
			+ ' | grep -v H'
			+ ' | grep -v score'
			+ ' | grep -v pdb'
			+ ' | grep -v descr'
			+ ' | grep -v total'
			+ ' | grep -v S_' 
	)
	scores = popen( command ).read().split('\n')[:-1]
	scores = [ score.split() for score in scores ]

	data.scores = scores
	data.score_labels = score_labels

	return data

##########################################################

def get_date():
	### Option 1: use datetime
	#from datetime import date
	#today = date.today()
	#current_date = "%d/%d/%d" % ( today.month, today.day, today.year )
	### Option 2: use popen( 'date' )
	current_date = popen( 'date' ).read().replace('\n','')
	return current_date

##########################################################

def get_target_names( target_files ):
	target_names = []
	for file_name in target_files:
		if '..' not in file_name:	file_name = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/input_files/' + file_name
		target_names = get_target_names_from_file( file_name, target_names )
	return target_names

###########################################################

def get_target_names_from_file( filename, target_names ):
	assert( exists( filename ) )
	fid = open( filename, 'r' )
	for line in fid.readlines():
		### Old target file format
		### Name 			Sequence 	Secstruct ...
		### target_name
		#if 'Name' in line: continue
		#cols = line.split()
		#if not len( cols ): continue
		#target_names.append( cols[0] )
		### New target file format
		### Name:		target_name
		cols = string.split( line.replace( '\n', '' ) )
		if not len( cols ): continue
		if cols[0] == 'Name:': target_names.append( cols[1] ) 
	fid.close()
	return target_names

###########################################################

def get_path_to_dir( dirnames ):
	for dirname in dirnames:
		pwd = popen( 'pwd' ).readline()[:-1].split( '/' )
		if dirname not in pwd: continue
		return string.join( pwd[:pwd.index( dirname )+1], '/' )

###########################################################

def show_times( inpaths, data, noutfiles, target_names, which_target ):
	times_list = get_times( inpaths, data, noutfiles, target_names, which_target )
	for k in xrange( len(target_names) ):
		times_found = False
		for n in xrange( len( inpaths ) ):
			if times_list[n][k].times_found():
				times_found = True
		if not times_found: continue

		print '\n %-6s%33s' % ( 'TARGET', target_names[k] )
		for n in xrange( len(inpaths) ):
			run_time = times_list[n][k]
			print ' Run %d                    %5.0f +/- %4.0f' % ( n+1, run_time.mean, run_time.stdev )
	print '\n'
	for n in xrange( len(inpaths) ):
		print ' Run %d: %s' % ( n+1, basename( inpaths[n] ) )
	return

###########################################################

def get_times( inpaths, data, noutfiles, target_names, which_target ):
	time_label = 'time'
	times_list = []
	for n in xrange( len( inpaths )):
		times = []
		for target in target_names:
			time_data = TimeData()
			time_data.inpath = inpaths[n]
			time_data.target = target
			try:
				time_idx = data[n][ target ].score_labels.index( time_label )
				time_data.times = np.array( [ score[time_idx] for score in data[n][ target ].scores ], dtype = 'float_' )
			except:
				time_data.times = [] 
			time_data.update_stats()
			times.append( time_data )
		times_list.append( times )
	return times_list

###########################################################

def get_figure_dimensions( noutfiles, landscape=False ):
	nplots = noutfiles
	assert( nplots )
	if landscape:
		if nplots < 3: 	  
			nrows = 1
		else:
			nrows = 3
	else:
		if nplots < 3: 	  
			nrows = nplots
		elif nplots < 12: 
			nrows = 4
		else:  		      
			nrows = 5
	ncols = np.ceil( nplots / float( nrows ) )

	if landscape:	
		figwidth = 11
		figheight = 8.5
	else:			
		figwidth = 8.5
		figheight = 11

	return ( nplots, nrows, ncols, figwidth, figheight ) 

###########################################################

def setup_pdf_page( base_inpaths, landscape=False, verbose=True ):
	pdfname = string.join(base_inpaths, '_vs_') 
	fullpdfname = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/Figures/' + pdfname 
	if landscape:	
		fullpdfname += '_landscape.pdf'
	else:			
		fullpdfname += '.pdf'
	print '\nMaking figure in: %s\n' % fullpdfname
	pp = PdfPages( fullpdfname )
	return ( pp, fullpdfname )

###########################################################

def jet( size ):
	#cmap = plt.get_cmap( 'jet' )
	cmap = plt.get_cmap( 'hot' )
	cnorm = colors.Normalize( vmin=0, vmax=size )
	scalar_map = cmx.ScalarMappable( norm=cnorm, cmap=cmap )
	colorcode = [ scalar_map.to_rgba( x ) for x in xrange( size ) ]
	return colorcode

###########################################################

def get_title( target ):
	names = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/python/plot_util/titles.txt'
	try:
		lines = open( names, 'r' ).readlines()
	except:
		return target
	for line in lines:
		if target in line:
			return line.replace('\n','').split(': ')[1]
	return target
