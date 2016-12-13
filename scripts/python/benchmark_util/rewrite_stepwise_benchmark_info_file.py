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

parser = argparse.ArgumentParser(description='Rewrite info text files with new format')
parser.add_argument('info_file', help='text file with motif definitions, in same directory as input_files/ (e.g., "input_files/favorites.txt")',default=None )
parser.add_argument('--tsv', help="convert to TSV file format", action="store_true")	
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

is_new_format = False


# make sure that info_file is specified correctly and exists
info_file = args.info_file
assert( len( info_file ) > 0 )
assert( '.txt' in info_file )
assert( exists( info_file )  )


# header for old format
header = [ 'Name','Sequence','Secstruct','Working_res','Native','Input_res','Extra_flags' ]
info_file_lines = open( info_file, 'r' ).readlines()
if string.split( info_file_lines[0].replace( '\n', '' ) ) != header: is_new_format = True


# read current info file
if not is_new_format:
	
	for info_file_line in info_file_lines[1:]:

		if info_file_lines[0] == '#':	continue
		if len( info_file_line ) < 5:	continue

		cols = string.split( info_file_line.replace( '\n', '' ) )

		name = cols[0]
		assert( name not in names )
		names.append( name )
		sequence   [ name ] = cols[1]
		secstruct  [ name ] = cols[2]
		working_res[ name ] = cols[3]
		native     [ name ] = cols[4]
		input_res  [ name ] = cols[5]
		extra_flags[ name ] = string.join( cols[6:] )

else:

	for info_file_line in info_file_lines:

		if info_file_lines[0] == '#':	continue
		if len( info_file_line ) < 5:	continue
		cols = string.split( info_file_line.replace( '\n', '' ) )

		if cols[0] == 'Name:' : 
			name = cols[1]
			assert( name not in names )
			names.append( name )

		if cols[0] == 'Sequence:'   : sequence   [ name ] = cols[1]
		if cols[0] == 'Secstruct:'  : secstruct  [ name ] = cols[1]
		if cols[0] == 'Working_res:': working_res[ name ] = cols[1]
		if cols[0] == 'Input_res:'  : input_res  [ name ] = cols[1]
		if cols[0] == 'Native:'     : native     [ name ] = cols[1]
		if cols[0] == 'Extra_flags:': extra_flags[ name ] = string.join( cols[1:] )


# check that each dictionary is the same size
assert( len( names ) == len( sequence ) == len( secstruct ) == len( working_res ) == len( input_res ) == len( native ))


# write new file
info_fid = open( info_file, 'w' )
if args.tsv:	
        lines = [ header ]
        for name in names:
                lines.append([
                        name,
                        sequence[ name ],
                        secstruct[ name ],
                        working_res[ name ],
                        native[ name ],
                        input_res[ name ],
                        extra_flags[ name ]
                ])  
        info_fid.write('\n'.join(map('\t'.join, lines)))
else:
	for name in names:
		info_fid.write( '%-15s%s\n' % ( 'Name:'       , name                ) )
		info_fid.write( '%-15s%s\n' % ( 'Sequence:'   , sequence   [ name ] ) )
		info_fid.write( '%-15s%s\n' % ( 'Secstruct:'  , secstruct  [ name ] ) )
		info_fid.write( '%-15s%s\n' % ( 'Working_res:', working_res[ name ] ) )
		info_fid.write( '%-15s%s\n' % ( 'Input_res:'  , input_res  [ name ] ) )
		info_fid.write( '%-15s%s\n' % ( 'Native:'     , native     [ name ] ) )
		info_fid.write( '%-15s%s\n' % ( 'Extra_flags:', extra_flags[ name ] ) )
		info_fid.write( '\n' )
info_fid.close()
