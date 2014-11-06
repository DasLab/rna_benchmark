#!/usr/bin/python

import argparse
from os.path import basename, dirname, exists
from os import popen, system, getcwd, chdir
import string
from read_pdb import read_pdb
from parse_tag import parse_tag
from get_surrounding_res import get_surrounding_res_tag

parser = argparse.ArgumentParser(description='Setup benchmark for input files for SWA')
parser.add_argument('info_file', help='text file with information, in same directory as input_files/ (e.g., "../favorites.txt")')
parser.add_argument('-r','--radius', help='Surrounding radius for input residues', default=None)
parser.add_argument('-s','--sequence',action='store_true')
parser.add_argument('--update_info_file',action='store_true')

args = parser.parse_args()


# read in benchmark information & create any missing input files.
lines = open( args.info_file ).readlines()
relpath = dirname( args.info_file )
if len( relpath ) > 0: relpath += '/'
header = lines[0]

# initialize directories
names = []
sequence = {}
secstruct = {}
working_res = {}
native = {}
input_res = {}
inpath = {}
extra_flags = {}
fasta = {}

if args.update_info_file:
	open( args.info_file, 'w' ).write( string.join( header.split(), '\t' ) + '\n' )

for line in lines[1:]:

	if line[0] == '#': continue
	cols = string.split( line.replace( '\n','' ) )
	name = cols[0]
	#print '\nSetting up input files for: ', name

	assert( not name in names ) # better be unique
	names.append( name )
	sequence [ name ]   = cols[1]
	secstruct[ name ]   = cols[2]
	working_res[ name ] = cols[3]
	native[ name ]      = cols[4]
	input_res[ name ]   = cols[5]
	inpath[ name ]      = relpath + 'input_files/' + basename( args.info_file.replace('.txt','' ) )
	extra_flags[ name ] = string.join( cols[6:] )

	
	if args.radius and input_res[ name ] == '-':
		try:
			surrounding_res_tag = get_surrounding_res_tag( inpath[ name ]+'/'+native[ name ], sample_res_list=working_res[ name ], radius=args.radius, verbose=False )
			input_res[ name ] = surrounding_res_tag.replace(' ',',')
			if input_res[ name ][0] == ',': input_res[ name ] = input_res[ name ][1:]
			#print 'Input_res: ',input_res[ name ]
		except:
			input_res[ name ] = '-'
	

	if args.sequence:
		if sequence[ name ] == '-':
			( coords, pdb_lines, sequence, chains, residues ) = read_pdb( inpath[ name ]+'/'+native[ name ] )
			if input_res[ name ] != '-':	( residues, chains ) = parse_tag( input_res[ name ]+' '+working_res[ name ] )
			seq = []
			for res_idx in xrange( len( residues ) ):
				seq.append( sequence[chains[res_idx]][residues[res_idx]])
			#if not len( seq ) > 40:
			sequence[ name ] = string.join( [res.replace( ' MG', '[MG]' ).replace(' ','').lower() for res in seq] ,'')
			#else:
			#	sequence[ name ] = '%s.fasta' % name
				### write fasta?
			#print 'Sequence: ',sequence[ name ]
   
   	print
	print 'Name:         ', name
	print 'Sequence:     ', sequence[ name ]
	print 'Secstruct:    ', secstruct[ name ] 
	print 'Working_res:  ', working_res[ name ] 
	print 'Native:       ', native[ name ] 
	print 'Input_res:    ', input_res[ name ] 
	print 'Extra_flags:  ', extra_flags[ name ] 


	info_line = [ name, sequence[ name ], secstruct[ name ], working_res[ name ], native[ name ], input_res[ name ], extra_flags[ name ] ]
	if args.update_info_file:
		open( args.info_file, 'a' ).write( string.join( info_line, '\t' ) + '\n' )



