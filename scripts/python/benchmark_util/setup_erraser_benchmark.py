#!/usr/bin/env python


import string
import argparse
from os.path import exists,basename,dirname,expandvars
from os import system, getcwd, chdir, popen
from setup_stepwise_benchmark_util import *
from sys import argv, exit
import subprocess

#####################################################################################################################

parser = argparse.ArgumentParser(description='Setup benchmark for ERRASER')
parser.add_argument("info_file",       help='text file with information, in same directory as input_files/ (e.g., "../*.erraser.txt")')
parser.add_argument("user_input_runs", nargs='*',help='specify particular cases to run (default: run all in info_file)' )
default_extra_flags_benchmark = 'extra_flags_benchmark.txt'
parser.add_argument('-extra_flags', default=default_extra_flags_benchmark, help='Filename of text file with extra_flags for all cases.')
parser.add_argument('-nhours', default='16', type=int, help='Number of hours to queue each job.')
parser.add_argument('-j','--njobs', default='1', type=int, help='Number of cores for each job.')
parser.add_argument('--path_to_rosetta', default='', help='Path to working copy of rosetta.')
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
parser.add_argument('--save_logs', help="save .out and .err logs for each job.", action="store_true")
args = parser.parse_args()

#####################################################################################################################


# get path to rosetta, required for now
ROSETTA=args.path_to_rosetta
if ( not len( ROSETTA ) ) or ( not exists( ROSETTA ) ):
    ROSETTA=expandvars( '$ROSETTA' )
    if ( not len( ROSETTA ) ) or ( not exists( ROSETTA ) ):
        print 'WARNING: $ROSETTA must be defined as the path to a working rosetta repository!!!'
        print 'Export this variable, by putting the following in your .bashrc or .zshrc:'
        print 'export ROSETTA=/path/to/rosetta/\n'
        exit(0)
assert( exists( ROSETTA ) )
ROSETTA_DB=ROSETTA+'/main/database/'
ERRASER_TOOLS=ROSETTA+'/tools/ERRASER/'


# get extra_flags_benchmark
if len (args.extra_flags) > 0:
    if exists( args.extra_flags ):
        extra_flags_benchmark = open( args.extra_flags ).readlines()
    else:
        extra_flags_benchmark = []
        print 'Did not find ', args.extra_flags, ' so not using any extra flags for the benchmark'
        assert ( args.extra_flags == default_extra_flags_benchmark )


# initialize directories
names = []
input_pdb = {}
map_file = {}
map_reso = {}
fixed_res = {}
extra_flags = {}


# make sure the info file is specified correctly and exists
info_file = args.info_file
assert( len( info_file ) > 0 )
assert( '.erraser.txt' in info_file )
assert( exists( info_file ) )


# define and check paths
inpath = info_file.replace('.erraser.txt', '' ) + '/'
assert( exists( inpath ) )


# read info_file
for info_file_line in open( info_file ).readlines():

    if info_file_line[0] == '#' : continue
    if len( info_file_line ) < 5: continue
    cols = string.split( info_file_line.replace( '\n','' ) )
    assert( len( cols ) >= 2 )

    if cols[0] == 'Name:':
        name = cols[1]
        if ( len( args.user_input_runs ) > 0 ) and ( name not in args.user_input_runs ): continue
        print 'Reading in info for: ', name
        assert( name not in names )
        names.append( name )
        continue

    if ( len( args.user_input_runs ) > 0 ) and ( name not in args.user_input_runs ): continue

    if   cols[0] == 'Input_pdb:'   :    input_pdb  [ name ] = cols[1]
    elif cols[0] == 'Map_file:'    :    map_file   [ name ] = cols[1]
    elif cols[0] == 'Map_reso:'    :    map_reso   [ name ] = cols[1]
    elif cols[0] == 'Fixed_res:'   :    fixed_res  [ name ] = cols[1].replace(':','').replace(',',' ') #erraser format A:33 --> A33 
    elif cols[0] == 'Extra_flags:' :    extra_flags[ name ] = string.join( cols[1:] )


# check that each dictionary is the same size
assert( len( names ) == len( input_pdb ) == len( map_file ) == len( map_reso ) == len( fixed_res ) )


# write qsubMINIs, READMEs and SUBMITs
qsub_file = 'qsubMINI'
hostname, hostname_err = subprocess.Popen(['hostname'], stdout=subprocess.PIPE).communicate()
if hostname.find( 'stampede' ) > -1: qsub_file = 'qsubMPI'
if hostname.find( 'sherlock' ) > -1: qsub_file = 'sbatchMINI'
fid_qsub = open( qsub_file, 'w' )

for name in names:

    dirname = name
    if not exists( dirname ): 
        system( 'mkdir '+dirname )

    # move all required files to the correct directory
    start_files = [ input_pdb[ name ], map_file[ name ] ]
    for start_file in start_files:
        system( 'cp %s/%s %s/ ' % ( inpath, start_file, dirname ) )

    # SETUP for ERRASER
    fid = open( '%s/README_ERRASER' % name, 'w' )
    fid.write( '%s/erraser.py' % ERRASER_TOOLS )
    fid.write( ' -pdb %s' % input_pdb[ name ] )
    fid.write( ' -map %s' % map_file[ name ] )
    if len(map_reso[ name ]) > 0 and map_reso[ name ] != '-':
        fid.write( ' -map_reso %s' % map_reso[ name ] )
    if len(fixed_res[ name ]) > 0 and fixed_res[ name ] != '-':
        fid.write( ' -fixed_res %s' % fixed_res[ name ] )

    # case-specific extra flags
    if ( len( extra_flags[name] ) > 0 ) and ( extra_flags[ name ] != '-' ) :
        for flag in extra_flags[name].split('-'):
            if not len( flag ): continue
            flag = flag.replace('true','True').replace('false','False')
            fid.write( ' -%s' % flag )

    for flag in extra_flags_benchmark:
        if ( '#' in flag ): continue
        flag = flag.replace('true','True').replace('false','False')
        fid.write( ' %s' % flag.replace('\n','') )

    fid.close()

    print '\nSetting up submission files for: ', name
    CWD = getcwd()
    chdir( name )
    rosetta_submit_cmd = 'rosetta_submit.py README_ERRASER OUT %d %d' % (args.njobs, args.nhours )
    if args.save_logs:
        rosetta_submit_cmd += ' -save_logs'
    system( rosetta_submit_cmd )
    chdir( CWD )
    fid_qsub.write( 'cd %s; source %s; cd %s\n' % ( name, qsub_file,  CWD ) )


fid_qsub.close()

