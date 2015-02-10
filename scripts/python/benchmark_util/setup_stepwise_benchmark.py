#!/usr/bin/env python


import string
import argparse
from os.path import exists,basename,dirname,expandvars
from os import system, getcwd, chdir, popen
from make_tag import *
from parse_options import get_resnum_chain
from parse_tag import parse_tag
from get_sequence import get_sequences
from rna_server_conversions import get_all_stems, join_sequence
from get_surrounding_res import get_surrounding_res_tag
from setup_stepwise_benchmark_util import *
from sys import argv, exit
import subprocess

#####################################################################################################################

parser = argparse.ArgumentParser(description='Setup benchmark for stepwise monte carlo')
parser.add_argument("info_file",       help='text file with information, in same directory as input_files/ (e.g., "../favorites.txt")')
parser.add_argument("user_input_runs", nargs='*',help='specify particular cases to run (default: run all in info_file)' )
default_extra_flags_benchmark = 'extra_flags_benchmark.txt'
parser.add_argument('-extra_flags', default=default_extra_flags_benchmark, help='Filename of text file with extra_flags for all cases.')
parser.add_argument('-nhours', default='16', type=int, help='Number of hours to queue each job.')
parser.add_argument('-j','--njobs', default='10', type=int, help='Number of cores for each job.')
parser.add_argument('--swa', action='store_true', help='Additional flag for setting up SWA runs.')
parser.add_argument('--extra_min_res_off', action='store_true', help='Additional flag for turning extra_min_res off.')
parser.add_argument('--save_times_off', action='store_true', help='Additional flag for turning save_times flag off.')
parser.add_argument('--path_to_rosetta', default='', help='Path to working copy of rosetta.')
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
parser.add_argument('-motif_mode_off', help="temporary hack for turning off hardcoded '-motif_mode' flag", action="store_true")
parser.add_argument('--save_logs', help="save .out and .err logs for each job.", action="store_true")
args = parser.parse_args()

#####################################################################################################################
motif_mode_off = ( args.motif_mode_off or args.swa )

if args.swa and args.njobs == 10:
    # set default njobs to 150 for SWA jobs
    njobs = 150
else:
    njobs = args.njobs


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
SWA_DAGMAN_TOOLS=ROSETTA+'/tools/SWA_RNA_python/SWA_dagman_python/'


# get extra_flags_benchmark
VDW_rep_screen_info_flag_found = False
cycles_flag_found = False
if len (args.extra_flags) > 0:
    if exists( args.extra_flags ):
        extra_flags_benchmark = open( args.extra_flags ).readlines()

        for flag in extra_flags_benchmark:
            if ( '-VDW_rep_screen_info' in flag ):
                VDW_rep_screen_info_flag_found = True
                continue
            if ( '-cycles' in flag ):
                cycles_flag_found = True
                continue

    else:
        extra_flags_benchmark = []
        print 'Did not find ', args.extra_flags, ' so not using any extra flags for the benchmark'
        assert ( args.extra_flags == default_extra_flags_benchmark )

# initialize directories
names = []
sequence = {}
secstruct = {}
working_res = {}
native = {}
input_res = {}
extra_flags = {}
fasta = {}
resnums = {}
chains  = {}
helix_files = {}
working_native = {}
input_pdbs = {}
terminal_res = {}
extra_min_res = {}
loop_res = {}
VDW_rep_screen_pdb = {}
VDW_rep_screen_info = {}
bps = ['au','ua','gc','cg','ug','gu']


# make sure the info file is specified correctly and exists
info_file = args.info_file
assert( len( info_file ) > 0 )
assert( '.txt' in info_file )
assert( exists( info_file ) )


# define and check paths
inpath = info_file.replace('.txt', '' ) + '/'
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

    if   cols[0] == 'Sequence:'    :    sequence   [ name ] = cols[1]
    elif cols[0] == 'Secstruct:'   :    secstruct  [ name ] = cols[1]
    elif cols[0] == 'Working_res:' :    working_res[ name ] = cols[1]
    elif cols[0] == 'Input_res:'   :    input_res  [ name ] = cols[1]
    elif cols[0] == 'Native:'      :    native     [ name ] = cols[1]
    elif cols[0] == 'Extra_flags:' :    extra_flags[ name ] = string.join( cols[1:] )


# check that each dictionary is the same size
assert( len( names ) == len( sequence ) == len( secstruct ) == len( working_res ) == len( input_res ) == len( native ))


# iterate over names
for name in names:

    sequences          = string.split( sequence[name], ',' )
    working_res_blocks = string.split( working_res[name], ',' )


    # store information on 'conventional' residue numbers and chains.
    resnums[ name ] = []
    chains[ name ] = []
    for working_res_block in working_res_blocks: get_resnum_chain( working_res_block, resnums[ name ], chains[ name ] )


    # working_native
    assert( native[ name ] != '-' ) # for now, require a native, since this is a benchmark.
    prefix = '%s/%s_' % ( inpath,name)
    working_native[ name ] = slice_out( inpath, prefix, native[ name ], string.join( working_res_blocks ) )
    assert( string.join(sequences,'') == string.join(get_sequences( working_native[name] )[0],'') )

    # create starting PDBs
    input_pdbs[ name ] = []
    input_resnums = []
    input_chains  = []
    if input_res[ name ] != '-':
        input_res_blocks = string.split( input_res[ name ], ';' )
        for m in range( len ( input_res_blocks ) ):
            prefix = '%s/%s_START%d_' % ( inpath,name,m+1)
            input_pdb = slice_out( inpath, prefix, native[ name ],input_res_blocks[m] )
            input_pdbs[ name ].append( input_pdb )
            get_resnum_chain( input_res_blocks[m], input_resnums, input_chains )


    input_resnum_fullmodel = map( lambda x: get_fullmodel_number(x,resnums[name],chains[name]), zip( input_resnums, input_chains ) )


    # create secstruct if not defined
    if secstruct[ name ] == '-':
        secstruct[ name ] = string.join( [ '.' * len( seq ) for seq in sequences ], ',' )


    # create any helices.
    helix_files[ name ] = []
    (sequence_joined, chainbreak_pos)           = join_sequence( sequence[name] )
    (secstruct_joined,chainbreak_pos_secstruct) = join_sequence( secstruct[name] )
    assert( chainbreak_pos == chainbreak_pos_secstruct )
    stems = get_all_stems( secstruct_joined, chainbreak_pos, sequence_joined  )

    for i in range( len( stems ) ):
        helix_file =  '%s/%s_HELIX%d.pdb' % (inpath,name,(i+1))

        stem = stems[i]
        helix_seq = ''; helix_resnum = [];
        for bp in stem:
            helix_seq    += sequence_joined[ bp[0] - 1 ]
            helix_resnum.append( bp[0] )
        helix_seq += ' '
        for bp in stem[::-1]:
            helix_seq    += sequence_joined[ bp[1] - 1 ]
            helix_resnum.append( bp[1] )

        already_in_input_res = False
        for m in helix_resnum:
            if m in input_resnum_fullmodel: already_in_input_res = True
        if already_in_input_res: continue

        helix_files[ name ].append( helix_file )
        input_resnum_fullmodel += helix_resnum

        if exists( helix_file ): continue
        command = 'rna_helix.py -seq %s  -o %s -resnum %s' % ( helix_seq, helix_file, \
            make_tag_with_conventional_numbering( helix_resnum, resnums[ name ], chains[ name ] ) )
        print command
        system( command )

    # following is now 'hard-coded' into Rosetta option '-motif_mode'
    # deprecate this python block in 2015 after testing -- rd2014
    L = len( sequence_joined )
    terminal_res[ name ] = []
    extra_min_res[ name ] = []
    for m in range( 1, L+1 ):
        if ( m not in input_resnum_fullmodel ): continue
        prev_moving = ( m - 1 not in input_resnum_fullmodel ) and ( m != 1 )
        next_moving = ( m + 1 not in input_resnum_fullmodel ) and ( m != L )
        right_before_chainbreak = ( m == L or m in chainbreak_pos )
        right_after_chainbreak  = ( m == 1 or m - 1 in chainbreak_pos )
        if ( ( right_after_chainbreak and not next_moving ) or \
             ( right_before_chainbreak and not prev_moving ) ):
            terminal_res[ name ].append( m )
        if ( ( prev_moving and not next_moving and not right_before_chainbreak ) or \
             ( next_moving and not prev_moving and not right_after_chainbreak ) ):
            extra_min_res[ name ].append( m )
    if not '-motif_mode\n' in extra_flags_benchmark and not motif_mode_off:
        extra_flags_benchmark.append( '-motif_mode\n' )

    # create fasta
    fasta[ name ] = '%s/%s.fasta' % (inpath,name)
    if not exists( fasta[ name ] ):
        fid = open( fasta[ name ], 'w' )
        assert( len( sequences ) == len( working_res_blocks ) )
        ### splitting up sequence in fasta may cause errors in SWA runs
        for n in range( len( sequences ) ): fid.write( '>%s %s\n%s\n' % (name,working_res_blocks[n],sequences[n]) )
        #fid.write( popen( 'pdb2fasta.py %s' % (  working_native[ name ] ) ).read() )
        #fid.write( '>%s %s\n%s\n' % ( name,string.join(working_res_blocks,' '),string.join(sequences,'') ) )
        fid.close()


    # get sample loop res
    loop_res[ name ] = {}
    if input_res[ name ] == '-':
        if args.swa:
            print "WARNING: input_res[ name ] == '-' "
        continue

    ( workres , workchains  ) = parse_tag( working_res[ name ], alpha_sort=True )
    ( inputres , inputchains  ) = parse_tag( input_res[ name ], alpha_sort=True )

    loopres_tag = []
    for ii in xrange( len( workres ) ):
        working_tag = workchains[ ii ] + ':' + str(workres[ ii ])
        is_input_tag = False
        for jj in xrange( len( inputres ) ):
            input_tag = inputchains[ jj ] + ':' + str(inputres[ jj ])
            if input_tag == working_tag:
                is_input_tag = True
        if is_input_tag: continue
        loopres_tag.append( working_tag )
    loopres_tag = string.join( loopres_tag, ',' )

    ( loopres , loopchains  ) = parse_tag( loopres_tag, alpha_sort=True )
    ( workres , workchains  ) = parse_tag( working_res[ name ], alpha_sort=True )

    loopres_conventional = [ str(workchains[idx])+':'+str(workres[idx]) for idx in xrange( len( workres ) ) if (workres[idx] in loopres and workchains[idx] == loopchains[loopres.index(workres[idx])]) ]
    loopres_conventional = string.join( [ str(x) for x in loopres_conventional ] ,' ')
    loop_res[ name ][ 'conventional' ] = loopres_conventional

    if args.swa:
        loopres_swa = [ idx+1 for idx in xrange( len( workres ) ) if (workres[idx] in loopres and workchains[idx] == loopchains[loopres.index(workres[idx])]) ]
        loopres_swa = string.join( [ str(x) for x in loopres_swa ] ,' ')
        loop_res[ name ][ 'swa' ]  = loopres_swa


    # get VDW_rep_screen_info, it will only be used if -VDW_rep_screen_info flag is set in extra_flags_benchmark
    periph_res_radius = 50.0
    if ( 'rrna' in name ) or ( 'rRNA' in name ):    periph_res_radius = 100.0

    prefix = '%s/%s_%d_ANGSTROM_GRID_' % ( inpath, name, periph_res_radius )
    VDW_rep_screen_pdb[ name ] = prefix + native[ name ]
    VDW_rep_screen_info[ name ] = '%s' % ( basename( VDW_rep_screen_pdb[name] ) )

    if VDW_rep_screen_info_flag_found:

        if not exists( VDW_rep_screen_pdb[ name ] ):

            loopres_list=string.split( loop_res[ name ][ 'conventional' ], ' ' )
            periph_res_tag = get_surrounding_res_tag( inpath+native[ name ], sample_res_list=loopres_list, radius=periph_res_radius, verbose=args.verbose )
            assert( len( periph_res_tag ) )
            slice_out( inpath, prefix, native[ name ], periph_res_tag )

            if args.verbose:
                print 'loopres_list for '+name+' = '+string.join(loopres_list)
                print 'periph_res for '+name+' = '+periph_res_tag



# write qsubMINIs, READMEs and SUBMITs
qsub_file = 'qsubMINI'
hostname = subprocess.check_output( 'hostname' )
if hostname.find( 'stampede' ) > 0: qsub_file = 'qsubMPI'
fid_qsub = open( qsub_file, 'w' )
for name in names:

    dirname = name
    if not exists( dirname ): system( 'mkdir '+dirname )

    # move all required files to the correct directory
    start_files = helix_files[ name ] + input_pdbs[ name ]
    infiles = start_files + [ fasta[name], working_native[ name ] ]
    if VDW_rep_screen_info_flag_found:  infiles.append( VDW_rep_screen_pdb[ name ] )
    for infile in infiles:  system( 'cp %s %s/ ' % ( infile, dirname ) )

    # SETUP for StepWise Assembly
    if args.swa:

        fid = open( '%s/README_SWA' % dirname, 'w' )
        fid.write( SWA_DAGMAN_TOOLS+'/SWA_DAG/setup_SWA_RNA_dag_job_files.py' )
        if len( start_files ) > 0 :
            fid.write( ' -s' )
            for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
        if len( native[ name ] ) > 0:
            fid.write( ' -native_pdb %s' % basename( working_native[name] ) )
        fid.write( ' -fasta %s.fasta' %  name )
        fid.write( ' -sample_res %s' % loop_res[ name ][ 'swa' ] )

        # case-specific extra flags
        if ( len( extra_flags[name] ) > 0 ) and ( extra_flags[ name ] != '-' ) :
            #fid.write( '%s\n' % extra_flags[name] )
            for flag in extra_flags[name].split('-'):
                if not len( flag ): continue
                flag = flag.replace('true','True').replace('false','False')
                fid.write( ' -%s' % flag )

        # extra flags for whole benchmark
        weights_file = ''
        for flag in extra_flags_benchmark:
            if ( '#' in flag ): continue
            flag = flag.replace('true','True').replace('false','False')
            if ( '-analytic_etable_evaluation' in flag ): continue ### SWM Specific
            if ( '-motif_mode' in flag ): continue ### SWM Specific
            if ( '-score:weights' in flag ):
                flag = flag.replace( '-score:weights', '-force_field_file' )
                weights_file = string.split( flag )[1]
            if ( '-score:rna_torsion_potential' in flag ):
                flag = flag.replace( '-score:rna_torsion_potential', '-rna_torsion_potential_folder' )
            if ( '-VDW_rep_screen_info True' in flag ):
                flag = flag.replace( 'True', VDW_rep_screen_info[ name ] ) #-VDW_rep_screen_info 1zih_RNA.pdb
                flag = flag + ' -apply_VDW_rep_delete_matching_res False'
            flag = ' '+flag.replace( '\n', '' )
            fid.write( flag )
        if len( weights_file ) > 0:
            if not exists( weights_file ):
                weights_file = ROSETTA_DB+'/scoring/weights/'+weights_file
            assert( exists( weights_file ) )
            system( 'cp ' + weights_file + ' ' + name )

        fid.close()

        print '\nSetting up submission files for: ', name
        CWD = getcwd()
        fid_submit = open( dirname+'/SUBMIT_SWA', 'w' )
        fid_submit.write( SWA_DAGMAN_TOOLS+'/dagman/submit_DAG_job.py' )
        fid_submit.write( ' -master_wall_time %d' % 72 ) #args.nhours )
        fid_submit.write( ' -master_memory_reserve 2048' )
        fid_submit.write( ' -num_slave_nodes %d' % njobs )
        fid_submit.write( ' -dagman_file rna_build.dag' )
        fid_submit.close()

        fid_qsub.write( 'cd %s; source ./README_SWA && source ./SUBMIT_SWA; cd %s\n' % ( dirname, CWD ) )

    # SETUP for StepWise Monte Carlo
    else:

        fid = open( '%s/README_SWM' % name, 'w' )
        fid.write( 'stepwise @flags -out:file:silent swm_rebuild.out\n' )
        fid.close()

        fid = open( '%s/flags' % name, 'w' )
        if len( start_files ) > 0 :
            fid.write( '-s' )
            for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
            fid.write( '\n' )
        if len( native[ name ] ) > 0:
            fid.write( '-native %s\n' % basename( working_native[name] ) )
        if len( terminal_res[ name ] ) > 0:
            fid.write( '-terminal_res %s  \n' % make_tag_with_conventional_numbering( terminal_res[ name ], resnums[ name ], chains[ name ] ) )
        if len( extra_min_res[ name ] ) > 0 and not args.extra_min_res_off: ### Turn extra_min_res off for SWM when comparing to SWA
            fid.write( '-extra_min_res %s \n' % make_tag_with_conventional_numbering( extra_min_res[ name ], resnums[ name ], chains[ name ] ) )
        #if ( len( input_pdbs[ name ] ) == 0 ):
        #    fid.write( '-superimpose_over_all\n' ) # RMSD over everything -- better test since helices are usually native
        fid.write( '-fasta %s.fasta\n' % name )
        if not cycles_flag_found:
            fid.write( '-cycles 200\n' )
        fid.write( '-nstruct 20\n' )
        #fid.write( '-intermolecular_frequency 0.0\n' )
        if not args.save_times_off:
            fid.write( '-save_times\n' )

        # case-specific extra flags
        if ( len( extra_flags[name] ) > 0 ) and ( extra_flags[ name ] != '-' ) :
            #fid.write( '%s\n' % extra_flags[name] )
            for flag in extra_flags[name].split('-'):
                if not len( flag ): continue
                flag = flag.replace('True','true').replace('False','false')
                fid.write( '-%s\n' % flag )

        # extra flags for whole benchmark
        weights_file = ''
        for flag in extra_flags_benchmark:
            if ( '#' in flag ): continue
            flag = flag.replace('True','true').replace('False','false')
            if ( '-single_stranded_loop_mode' in flag ): continue ### SWA Specific
            if ( '-score:weights' in flag ): weights_file = string.split( flag )[1]
            if ( '-VDW_rep_screen_info true' in flag ):
                flag = flag.replace( 'true', basename( VDW_rep_screen_info[ name ] ) )#-VDW_rep_screen_info 1zih_RNA.pdb
            fid.write( flag )

        if len( weights_file ) > 0:
            if not exists( weights_file ):
                weights_file = ROSETTA_DB+'/scoring/weights/'+weights_file
            assert( exists( weights_file ) )
            system( 'cp ' + weights_file + ' ' + name )

        fid.close()

        print '\nSetting up submission files for: ', name
        CWD = getcwd()
        chdir( name )

        rosetta_submit_cmd = 'rosetta_submit.py README_SWM SWM %d %d' % (njobs, args.nhours ) 
        if args.save_logs:
            rosetta_submit_cmd += ' -save_logs'
        system( rosetta_submit_cmd )

        chdir( CWD )

        fid_qsub.write( 'cd %s; source %s; cd %s\n' % ( name, qsub_file,  CWD ) )


fid_qsub.close()

