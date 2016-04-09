#!/usr/bin/python

from sys import argv
import argparse
from os import system
from os.path import exists
from make_tag import make_tag_with_dashes
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='generate a bunch of helices')
parser.add_argument('-weights', default='rna/rna_helix', help='Weights file defining score function')
parser.add_argument('-put_intra_into_total', action='store_true',default=False, help='calculate intra-res terms and include in totals')
args = parser.parse_args()

# create every 2bp and 3bp helix, and use to estimate 'turner rules'.
rna_seq = 'acgu'
rc = { 'a':'u', 'c':'g', 'g':'c', 'u':'a' }

def reverse_complement( seq ):
    seq_rc = ''
    for m in seq[::-1]:
        seq_rc += rc[ m ]
    return seq_rc

outdir = 'pdb'
if not exists( outdir ): system( 'mkdir ' + outdir )

commands = []
seq_sets = []

bps = ['au','ua','cg','gc','gu','ug'];

# singlets
for bp in bps:
    seq1 = bp[0]
    seq2 = bp[1]
    seq_sets.append( [seq1, seq2 ] )

# singlets with bulge
for m in rna_seq:
    seq1 = 'gc'+m
    seq2 = '0gc'
    seq_sets.append( [seq1, seq2 ] )

for m in rna_seq:
    seq1 = 'cg'+m
    seq2 = '0cg'
    seq_sets.append( [seq1, seq2 ] )

for m in rna_seq:
    seq1 = m+'cg'
    seq2 = 'cg0'
    seq_sets.append( [seq1, seq2 ] )

for m in rna_seq:
    seq1 = m+'gc'
    seq2 = 'gc0'
    seq_sets.append( [seq1, seq2 ] )

# doublets
for bp1 in bps:
    for bp2 in bps:
        seq1 = bp1[0] + bp2[0]
        seq2 = bp2[1] + bp1[1]

        seq_sets.append( [seq1, seq2 ] )


seq_sets.append( [ 'ggg', 'ccc' ] )
seq_sets.append( [ 'gaa', 'uuc' ] )
seq_sets.append( [ 'ccg', 'cgg' ] )
seq_sets.append( [ 'cca', 'ugg' ] )
seq_sets.append( [ 'aca', 'ugu' ] )

seq_sets.append( [ 'ggc', 'guc' ] )
seq_sets.append( [ 'cuc', 'ggg' ] )
# specially stable...
seq_sets.append( [ 'gguc', 'gguc' ] )
# new
seq_sets.append( [ 'cgcg', 'cgcg' ] )
seq_sets.append( [ 'cccc', 'gggg' ] )
seq_sets.append( [ 'acgu', 'acgu' ] )

for seq_set in seq_sets:

    seq1 = seq_set[0]
    seq2 = seq_set[1]

    chains = []
    resnum = []
    count = 0
    for i in range( len( seq1 ) ):
        if seq1[i] == '0': continue
        count += 1
        chains.append( 'A' )
        resnum.append( count )
    for i in range( len( seq2 ) ):
        if seq2[i] == '0': continue
        count += 1
        chains.append( 'B' )
        resnum.append( count )

    outfile = '%s/%s_%s_helix.pdb' % ( outdir, seq1, seq2 )
    if not exists( outfile ):
        command = 'rna_helix.py -seq %s %s -finish_weights %s -o %s -resnum %s ' % \
            (seq1, seq2, args.weights, outfile, make_tag_with_dashes( resnum, chains ) )
        if args.put_intra_into_total: command += '-put_intra_into_total '
        print command
        # system command
        commands.append( command )

print commands
pool = Pool(processes=4)              # start 4 worker processes
pool.map( system, commands )

# create every 2bp helix with 3' dangling ends.

# create every 2bp helix with 5' dangling ends.


