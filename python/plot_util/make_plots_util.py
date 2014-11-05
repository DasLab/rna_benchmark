#!/usr/bin/python

##########################################################

from os import system, popen
import string
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

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
		if '..' not in file_name:	file_name = get_path_to_dir('stepwise_benchmark') + '/' + file_name
		target_names = get_target_names_from_file( file_name, target_names )
	return target_names

###########################################################

def get_target_names_from_file( filename, target_names ):
	fid = open( filename, 'r' )
	for line in fid.readlines():
		if 'Name' in line: continue
		cols = line.split()
		if not len( cols ): continue
		target_names.append( cols[0] ) 
	fid.close()
	return target_names

###########################################################

def get_path_to_dir( dirname ):
	pwd = popen( 'pwd' ).readline()[:-1].split( '/' )
	return string.join( pwd[:pwd.index( dirname )+1], '/' )

###########################################################

def jet( size ):
	cmap = plt.get_cmap( 'hot' )
	cnorm = colors.Normalize( vmin=0, vmax=size )
	scalar_map = cmx.ScalarMappable( norm=cnorm, cmap=cmap )
	colorcode = [ scalar_map.to_rgba( x ) for x in xrange( size ) ]
	return colorcode