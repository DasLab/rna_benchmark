#!/usr/bin/python


import argparse
from os.path import basename, dirname, exists
from os import popen, system, getcwd, chdir
import string
from read_pdb import read_pdb
from parse_tag import parse_tag
from make_tag import make_tag_with_dashes, make_tag_with_dashes_and_commas
from get_surrounding_res import get_surrounding_res_tag


#####################################################################################################################

parser = argparse.ArgumentParser(description='Setup motif info text files to be read by setup_stepwise_benchmark.py')
parser.add_argument('loops_file', help='text file with loop motif definitions, in same directory as input_files/ (e.g., "input_files/favorites.loops")',default=None )
parser.add_argument('-radius', help='Surrounding radius for input residues', default=None)
parser.add_argument('--overwrite', help="overwrite existing info file", action="store_true")	
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")	
args = parser.parse_args()

#####################################################################################################################


# initialize directories
names = []
sequence = {}
secstruct = {}
working_res = {}
native = {}
input_res = {}
extra_flags = {}
fasta = {}
loop_res = {}


# make sure that loops_file is specified correctly and exists
loops_file = args.loops_file
assert( len( loops_file ) > 0 )
assert( '.loops' in loops_file )
assert( exists( loops_file )  )


# create info_file from loops_file and make sure it does not already exist
info_file = loops_file.replace( '.loops', '.txt' )
assert( '.txt' in info_file )
if not args.overwrite:	assert( not exists( info_file ) )
info_fid = open( info_file, 'w' )
info_fid.close()
assert( exists( info_file ) )


# define paths and check 
inpath = info_file.replace('.txt', '' ) + '/'
assert( exists( inpath ) )


# read in loop motifs information & get any missing information.
if args.verbose:	print '\nReading file: %s' % loops_file
for loops_file_line in open( loops_file ).readlines():
	

	if loops_file_line[0] == '#' : continue
	if len( loops_file_line ) < 5: 
		print
		continue 


	cols = string.split( loops_file_line.replace( '\n', '' ) )
	assert( len( cols ) == 2 )
	

	if cols[0] == 'Name:':
		name = cols[1]
		assert( not name in names )
		names.append( name )
		if args.verbose:	print '%-15s%s' % ( 'Name:', name )
		continue
		

	if cols[0] == 'Loop_res:':
		loop_res[ name ] = cols[1]
		if args.verbose:	print '%-15s%s' % ( 'Loop_res:', loop_res[ name ] )
		continue

	
	if cols[0] == 'Native:':
		native[ name ] = cols[1]
		assert( exists( inpath+native[ name ] ) )
		if args.verbose:	print '%-15s%s' % ( 'Native:', native[ name ] )
		continue
	

# check that each dictionary is the same size
assert( len( names ) == len( loop_res ) == len( native ) )


# write info file
if args.verbose:	print '\nWriting file: %s' % info_file
for name in names:

	
	# initialize values for key = name
	sequence   [ name ] = '-'
	secstruct  [ name ] = '-'
	working_res[ name ] = '-'
	input_res  [ name ] = '-'
	extra_flags[ name ] = '-'


	# get input res
	input_res[ name ] = get_surrounding_res_tag( inpath+native[ name ], sample_res_list=loop_res[ name ], radius=args.radius, csv=True )
	assert( len( input_res[ name ] ) )


	# get working res
	( inputres, inputchains ) = parse_tag( input_res[ name ] )
	( loopres , loopchains  ) = parse_tag( loop_res [ name ] )
	sorted_residues = zip(  ( inputchains + loopchains ), ( inputres + loopres ) )#.sort()
	sorted_residues.sort()
	[ working_chains, working_residues ] = [ list(x) for x in zip( *sorted_residues ) ]
	working_res[ name ] = make_tag_with_dashes_and_commas( working_residues, working_chains )
	
	### NEED TO REORDER WORKING_RES somehow

	assert( len( working_res[ name ] ) )


	# get sequence 
	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( inpath+native[ name ] )
	working_sequences = []
	working_res_blocks = working_res[ name ].split(',')
	for working_res_block in working_res_blocks:
		working_sequence = ''
		( working_residues, working_chains ) = parse_tag( working_res_block )
		for res_idx in xrange( len( working_residues ) ):
			xres   = working_residues[ res_idx ]
			xchain = working_chains  [ res_idx ]
			xseq   = sequence[ xchain ][ xres ].replace( ' ', '' ).lower()
			assert( xseq in [ 'a', 'u', 'c', 'g' ] )
			working_sequence += xseq
		assert( len( working_sequence ) )
		working_sequences.append( working_sequence )
	sequence[ name ] = string.join( working_sequences, ',' )
		 

	# print output to terminal if verbose flag on
	if args.verbose:
		print '%-15s%s' % ( 'Name:'       , name                )
		print '%-15s%s' % ( 'Sequence:'   , sequence   [ name ] )
		print '%-15s%s' % ( 'Secstruct:'  , secstruct  [ name ] )
		print '%-15s%s' % ( 'Working_res:', working_res[ name ] )
		print '%-15s%s' % ( 'Input_res:'  , input_res  [ name ] )
		print '%-15s%s' % ( 'Native:'     , native     [ name ] )
		print '%-15s%s' % ( 'Extra_flags:', extra_flags[ name ] )
		print
	
	# setup/write info_file
	info_fid = open( info_file, 'a' )
	info_fid.write( '%-15s%s\n' % ( 'Name:'       , name                ) )
	info_fid.write( '%-15s%s\n' % ( 'Sequence:'   , sequence   [ name ] ) )
	info_fid.write( '%-15s%s\n' % ( 'Secstruct:'  , secstruct  [ name ] ) )
	info_fid.write( '%-15s%s\n' % ( 'Working_res:', working_res[ name ] ) )
	info_fid.write( '%-15s%s\n' % ( 'Input_res:'  , input_res  [ name ] ) )
	info_fid.write( '%-15s%s\n' % ( 'Native:'     , native     [ name ] ) )
	info_fid.write( '%-15s%s\n' % ( 'Extra_flags:', extra_flags[ name ] ) )
 	info_fid.write( '\n' )
	info_fid.close()


