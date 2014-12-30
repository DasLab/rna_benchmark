#!/usr/bin/python

##########################################################

from os import system, popen
from os.path import exists, dirname, basename, abspath
import string
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

	command = 	'grep '+keyword+' '+file 								 	\
				+ ' | awk \'{print $'+string.join( scorenums, ',$' )+'}\''  \
				+ ' | grep -v inp' 											\
				+ ' | grep -v R' 											\
				+ ' | grep -v H' 											\
				+ ' | grep -v score' 										\
				+ ' | grep -v pdb' 											\
				+ ' | grep -v descr'										\
				+ ' | grep -v total'										\
				+ ' | grep -v S_' 											\

	scores = popen( command ).read().split('\n')[:-1]
	scores = [ score.split() for score in scores ]

	data.scores = scores
	data.score_labels = score_labels

	return data

##########################################################

def get_target_names( target_files ):
	target_names = []
	for file_name in target_files:
		if '..' not in file_name:	file_name = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/input_files/' + file_name
		target_names = get_target_names_from_file( file_name, target_names )
	return target_names

###########################################################

def get_target_names_from_file( filename, target_names ):
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

	time_name = 'time'
	times_list = []

	for n in xrange( len( inpaths ) ):
		times = []
		for target in target_names:
			try:
				time_idx = data[n][ target ].score_labels.index( time_name )
				times.append( np.array( [ score[time_idx] for score in data[n][ target ].scores ], dtype = 'float_' ) )
			except:
				times.append( [] )
                                #times.append( np.array( [ 0.0 for score in xrange( noutfiles ) ], dtype = 'float_' ) )
		times_list.append( times )


	for k in xrange( len(target_names) ):
                times_found = False
		for n in xrange( len( inpaths ) ):
                        if ( len( times_list[n][k] ) > 0.0 ): times_found = True
                if not times_found: continue

		print '\n %-6s%33s' % ( 'TARGET', target_names[ k ] ) #which_target[0][k] ] )
		for n in xrange( len( inpaths ) ):
			mean_time = np.mean( times_list[n][k] )
			std_time = np.std( times_list[n][k] )
			print ' Run %d                    %5.0f +/- %4.0f' % ( n, mean_time, std_time )

	print '\n'
	for n in xrange( len( inpaths ) ):
		print ' Run %d: %s' % ( n, basename( inpaths[n] ) )

	return

###########################################################

def jet( size ):
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
