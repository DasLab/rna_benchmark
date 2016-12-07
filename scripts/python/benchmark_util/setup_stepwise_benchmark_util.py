#!/usr/bin/python

#####################################################################################################################

import string
from os.path import exists
from os import system
from parse_options import get_resnum_chain
from make_tag import make_tag_with_dashes
from get_sequence import get_sequences
#####################################################################################################################

def flatten( l ):
    new_l = []
    for e in l: new_l.extend( e )
    return new_l

# helper function for PDB processing
def slice_out( inpath_dir, prefix, pdb, res_string, excise=False, check_sequence=False ):
    starting_native = inpath_dir+'/'+pdb
    assert( exists( starting_native ) )
    slice_pdb = prefix + pdb
    if not exists( slice_pdb ):
        if excise:  command = 'pdbslice.py %s -excise %s %s ' % ( starting_native, res_string, prefix )
        else:       command = 'pdbslice.py %s -subset %s %s ' % ( starting_native, res_string, prefix )
        system( command )
    assert( exists( slice_pdb ) )
    if check_sequence:
        target_resnums = []
        target_chains = []
        for col in res_string.split(' '):
            get_resnum_chain( col, target_resnums, target_chains )
        ( sequences, all_chains, all_resnums ) = get_sequences( slice_pdb )
        assert( sorted(flatten(all_resnums)) == sorted(target_resnums) )
        assert( sorted(flatten(all_chains)) == sorted(target_chains) )

    return slice_pdb

#####################################################################################################################

def get_fullmodel_number( reschain, resnums, chains):
    for m in range( len( resnums ) ):
        if ( resnums[m] == reschain[0] ) and ( reschain[1] == '' or chains[m] == reschain[1] ): return m+1
    return 0

#####################################################################################################################

def get_align_res( screen_pdb, working_pdb, working_fixed_res ):
    from read_pdb import read_pdb
    screen_align_res = []
    working_align_res = []
    for pdb in [ screen_pdb, working_pdb ]: assert( exists( pdb ) )
    #( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdb )
    screen_pdb_info = read_pdb( screen_pdb )
    working_pdb_info = read_pdb( working_pdb )
    for n in xrange( len( screen_pdb_info[4] ) ):
        screen_chain = screen_pdb_info[3][n]
        screen_res = screen_pdb_info[4][n]
        #working_reschain = pdb2pose( working_pdb_info[3], working_pdb_info[4], screen_chain, screen_res )
        working_reschain = ( 0, '' )
        for m in xrange( len( working_pdb_info[3] ) ):
            if ( working_pdb_info[4][m] == screen_res ) and ( working_pdb_info[3][m] == screen_chain ):
                working_reschain = ( working_pdb_info[4][m], working_pdb_info[3][m] )
        if ( working_reschain[0] > 0 ):# and ( working_reschain[0] in working_fixed_res ):
            screen_align_res.append( n )
            working_align_res.append( get_fullmodel_number(working_reschain,working_pdb_info[4],working_pdb_info[3]) )
    if len( screen_align_res ): screen_align_res_tag = make_tag_with_dashes( screen_align_res )
    else:   screen_align_res_tag = '0-0'
    if len( working_align_res ):    working_align_res_tag = make_tag_with_dashes( working_align_res )
    else:   working_align_res_tag = '0-0'
    return ( screen_align_res_tag, working_align_res_tag )
