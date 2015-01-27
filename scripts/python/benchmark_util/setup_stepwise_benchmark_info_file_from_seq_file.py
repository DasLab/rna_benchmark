#!/usr/bin/python


import argparse
from os.path import basename, dirname, exists
from os import popen, system, getcwd, chdir
import string
from read_pdb import read_pdb
from parse_tag import parse_tag
from make_tag import make_tag_with_dashes, make_tag_with_dashes_and_commas
from get_surrounding_res import get_surrounding_res_tag
import subprocess


#####################################################################################################################

parser = argparse.ArgumentParser(description='Setup motif info text files to be read by setup_stepwise_benchmark.py')
parser.add_argument('seq_file', help='text file with names and sequences, in same directory as input_files/ (e.g., "input_files/favorites.loops")',default=None )
parser.add_argument('-native_template', help='initial native pdb for threading sequence and creating unique native pdbs', default=None)
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
init_sequence = {}


# make sure that seq_file is specified correctly and exists
seq_file = args.seq_file
assert( len( seq_file ) > 0 )
assert( '.seq' in seq_file )
assert( exists( seq_file )  )


# create info_file from seq_file and make sure it does not already exist
info_file = seq_file.replace( '.seq', '.txt' )
assert( '.txt' in info_file )
assert( not exists( info_file ) or args.overwrite )
info_fid = open( info_file, 'w' )
info_fid.close()
assert( exists( info_file ) )


# define paths and check 
inpath = info_file.replace('.txt', '' ) + '/'
assert( exists( inpath ) )


# read in loop motifs information & get any missing information.
if args.verbose:	print '\nReading file: %s' % seq_file
for seq_file_line in open( seq_file ).readlines():

    if seq_file_line[0] == '#' : continue
    seq_file_line = seq_file_line.split('#')[0] # remove comments

    cols = string.split( seq_file_line.replace( '\n','' ) )
    assert( len( cols ) == 2 )

    name = cols[0]
    assert( name not in names )
    
    names.append( name )
    init_sequence[ name ] = cols[1]

    if args.verbose:
    	print name
    	print init_sequence[ name ]
    	print 

# check that each dictionary is the same size
assert( len( names ) == len( init_sequence ) )


# write info file
if args.verbose:	print '\nWriting file: %s' % info_file
for name in names:

	
	# initialize values for key = name
	sequence   [ name ] = '-'
	secstruct  [ name ] = '-'
	working_res[ name ] = '-'
	input_res  [ name ] = '-'
	native     [ name ] = '-'
	extra_flags[ name ] = '-'



	# hardcoding sequence reading for tandemGA
	## get start stop from init_native??
	start1, stop1 = 31, 38
	start2, stop2 = 45, 52
	sliced_sequences = [ init_sequence[ name ][start1:stop1+1], init_sequence[ name ][start2:stop2+1] ]
	sequence[ name ] = string.join( sliced_sequences , ',').lower()

	# hardcoding secstruct for tandemGA
	secstruct[ name ] = '(((..(((,)))..)))'

	# get native using rna_thread
	native_template = args.native_template
	native[ name ] = basename( native_template.replace( '.pdb', '_%s.pdb' % name ) )
	if not exists( native[ name ] ):
		rna_thread_cmdline = ['rna_thread', '-s', native_template, '-seq', sequence[ name ].replace(',','') , '-o', native[ name ] ] 
		#print 'running: ', string.join(rna_thread_cmdline, ' ')
		out, err = subprocess.Popen( rna_thread_cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()

	# get working res
	native_pdb_info = read_pdb( native[ name ] ) # ( coords, pdb_lines, sequence, chains, residues )
	working_chains = native_pdb_info[3]
	working_residues = native_pdb_info[4]
	working_res[ name ] = make_tag_with_dashes_and_commas( working_residues, working_chains )

	
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


