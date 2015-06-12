#!/usr/bin/env python


import string
import argparse
import os
import sys
from os.path import exists,basename,dirname,expandvars
from make_tag import *
from parse_options import get_resnum_chain
from parse_tag import parse_tag
from get_sequence import get_sequences
from rna_server_conversions import get_all_stems, join_sequence
from get_surrounding_res import get_surrounding_res_tag
from utility import helpers, info_handlers, file_handlers
from setup_stepwise_benchmark_util import *
from sys import argv, exit
import subprocess

#####################################################################################################################

parser = argparse.ArgumentParser(description='Setup benchmark for stepwise monte carlo')
parser.add_argument("info_file",       help='text file with information, in same directory as input_files/ (e.g., "../favorites.txt")')
parser.add_argument("user_input_runs", nargs='*',help='specify particular cases to run (default: run all in info_file)' )
parser.add_argument('-extra_flags', default='extra_flags_benchmark.txt', help='Filename of text file with extra_flags for all cases.')
parser.add_argument('-nhours', default='16', type=int, help='Number of hours to queue each job.')
parser.add_argument('-j','--njobs', default=None, type=int, help='Number of cores for each job.')
parser.add_argument('--extra_min_res_off', action='store_true', help='Additional flag for turning extra_min_res off.')
parser.add_argument('--save_times_off', action='store_true', help='Additional flag for turning save_times flag off.')
parser.add_argument('-motif_mode_off', help="temporary hack for turning off hardcoded '-motif_mode' flag", action="store_true")
parser.add_argument('--save_logs', help="save .out and .err logs for each job.", action="store_true")
parser.add_argument('--rosetta', default='', help='Path to working copy of rosetta.')
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
parser.add_argument('--swa', action='store_true', help='Additional flag for setting up SWA runs.')
parser.add_argument('--design', help="design all working residues not in input structures", action="store_true")
args = parser.parse_args()

#####################################################################################################################
motif_mode_off = ( args.motif_mode_off or args.swa )
njobs = 150 if args.swa else 10 
njobs = njobs if args.njobs is None else args.njobs 

# get path to rosetta, required for now
helpers.init_environ(args)

# replace python/c++ syntax accordingly
replacements = { 'True' : 'true', 'False' : 'false' }
if args.swa:
    replacements = { 'true' : 'True', 'false' : 'False' }


# parse extra_flags_benchmark
extra_flags_benchmark = {}
input_res_benchmark = None
extra_min_res_benchmark = None
if exists(args.extra_flags):
    with open(args.extra_flags, 'r') as fid:
        lines = filter(None, [l.split('#')[0].strip() for l in fid])
    extra_flags_benchmark = parse_flags(lines, replacements)
    if '-input_res' in extra_flags_benchmark:
        input_res_benchmark = extra_flags_benchmark.pop('-input_res')
    if '-extra_min_res' in extra_flags_benchmark:
        extra_min_res_benchmark = extra_flags_benchmark.pop('-extra_min_res')
else:
    print args.extra_flags,"doesn't exist, not using extra flags for benchmark"


# initialize directories
targets = []
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
info_fid = file_handlers.TargetDefinitionsFile()
info_fid.load(open(info_file))
assert( info_fid.validate() )
target_definitions = info_fid.target_definitions
if len(args.user_input_runs):
    target_definitions = [td for td in target_definitions if td.name not in args.user_input_runs]
targets = [info_handlers.Target(td) for td in target_definitions]


# misc dicts
full_model_info = {}

# iterate over names
for target in targets:

    full_model_info[ target.name ] = FullModelInfo( target.name )
    
    sequences          = string.split( target.sequence, ',' )
    working_res_blocks = string.split( target.working_res, ',' )

    # store information on 'conventional' residue numbers and chains.
    resnums = []
    chains = []
    for working_res_block in working_res_blocks:
        get_resnum_chain( working_res_block, resnums, chains )
    full_model_info[ target.name ].set_resnums( resnums )
    full_model_info[ target.name ].set_chains( chains )
    target.chains = chains
    target.resnums = resnums

    # working_native
    assert( target.native != '-' ) # for now, require a native, since this is a benchmark.
    prefix = '%s/%s_' % ( inpath, target.name)
    target.working_native = slice_out( inpath, prefix, target.native, string.join( working_res_blocks ) )
    assert( string.join(sequences,'') == string.join(get_sequences( target.working_native )[0],'') )

    # create starting PDBs
    target.input_pdbs = []
    input_res_blocks = [] 
    if target.input_res != '-':
        input_res_blocks += target.input_res.split(';')
    if input_res_benchmark:
        input_res_blocks += input_res_benchmark.split(';')
    for m,input_res_block in enumerate(input_res_blocks):
        prefix = '%s/%s_START%d_' % ( inpath,target.name,m+1)
        input_pdb = slice_out( inpath, prefix, target.native,input_res_block )
        target.input_pdbs.append( input_pdb )
    input_resnum_fullmodel = full_model_info[ target.name ].conventional_tag_to_full( 
        input_res_blocks
    )
    input_resnum_fullmodel.sort()

    #input_resnum_fullmodel = map( lambda x: get_fullmodel_number(x,resnums[name],chains[name]), zip( input_resnums, input_chains ) )


    # create secstruct if not defined
    if target.secstruct == '-':
        target.secstruct = string.join( [ '.' * len( seq ) for seq in sequences ], ',' )


    # create any helices.
    target.helix_files = []
    (sequence_joined, chainbreak_pos)           = join_sequence( target.sequence )
    (secstruct_joined,chainbreak_pos_secstruct) = join_sequence( target.secstruct )
    assert( chainbreak_pos == chainbreak_pos_secstruct )
    stems = get_all_stems( secstruct_joined, chainbreak_pos, sequence_joined  )

    for i in range( len( stems ) ):
        helix_file =  '%s/%s_HELIX%d.pdb' % (inpath,target.name,(i+1))

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

        target.helix_files.append( helix_file )
        input_resnum_fullmodel += helix_resnum
        input_resnum_fullmodel.sort()

        if exists( helix_file ): continue
        command = 'rna_helix.py -seq %s  -o %s -resnum %s' % ( helix_seq, helix_file, \
            make_tag_with_conventional_numbering( helix_resnum, resnums, chains) )
        print command
        os.system( command )
    

    # following is now 'hard-coded' into Rosetta option '-motif_mode'
    # deprecate this python block in 2015 after testing -- rd2014
    L = len( sequence_joined )
    target.terminal_res = []
    target.extra_min_res = []
    for m in range( 1, L+1 ):
        if ( m not in input_resnum_fullmodel ): continue
        prev_moving = ( m - 1 not in input_resnum_fullmodel ) and ( m != 1 )
        next_moving = ( m + 1 not in input_resnum_fullmodel ) and ( m != L )
        right_before_chainbreak = ( m == L or m in chainbreak_pos )
        right_after_chainbreak  = ( m == 1 or m - 1 in chainbreak_pos )
        if ( ( right_after_chainbreak and not next_moving ) or \
             ( right_before_chainbreak and not prev_moving ) ):
            target.terminal_res.append( m )
        if ( ( prev_moving and not next_moving and not right_before_chainbreak ) or \
             ( next_moving and not prev_moving and not right_after_chainbreak ) ):
            target.extra_min_res.append( m )
    if '-extra_min_res' in target.extra_flags:
        extra_res_full = full_model_info[ target.name ].conventional_tag_to_full(
            target.extra_flags.pop('-extra_min_res' )
        )
        extra_res_full = filter(input_resnum_fullmodel.count, extra_res_full)
        target.extra_min_res += extra_res_full
        target.extra_min_res = sorted(list(set(target.extra_min_res))) 
        motif_mode_off = True
    if extra_min_res_benchmark:
        extra_res_full = full_model_info[ target.name ].conventional_tag_to_full(
            extra_min_res_benchmark
        )
        extra_res_full = filter(input_resnum_fullmodel.count, extra_res_full)
        target.extra_min_res += extra_res_full
        target.extra_min_res = sorted(list(set(target.extra_min_res))) 
        motif_mode_off = True
    if not motif_mode_off and '-motif_mode' not in extra_flags_benchmark:
        extra_flags_benchmark['-motif_mode'] = ''

    # create fasta
    target.fasta = '%s/%s.fasta' % (inpath,target.name)
    if args.design:
        target.fasta = target.fasta.replace('.fasta', '_design.fasta')
        sequence_joined, chainbreak_pos = join_sequence( target.sequence )
        design_sequences = ['']
        for res_idx, res in enumerate(sequence_joined, start=1):
            if not res_idx in input_resnum_fullmodel:
                res = 'n'
            if res_idx - 1 in chainbreak_pos:
                design_sequences.append('')
            design_sequences[-1] += res
        sequences = design_sequences
        target.sequence = ','.join(design_sequences)
        print input_resnum_fullmodel
        print sequences
    if args.swa:
        target.fasta = target.fasta.replace('.fasta', '_SWA.fasta')
    if not exists( target.fasta ):
        fid = open( target.fasta, 'w' )
        assert( len( sequences ) == len( working_res_blocks ) )
        if args.swa:
            fid.write( '>%s %s\n%s\n' % ( target.name,string.join(working_res_blocks,' '),string.join(sequences,'') ) )
        else:
            ### splitting up sequence in fasta may cause errors in SWA runs
            for n in range( len( sequences ) ): fid.write( '>%s %s\n%s\n' % (target.name,working_res_blocks[n],sequences[n]) )
        #fid.write( os.popen( 'pdb2fasta.py %s' % (  working_native[ name ] ) ).read() )
        fid.close()

    # get align_pdb
    if '-align_pdb' in extra_flags_benchmark:
        prefix = '%s/%s_ALIGN_' % (inpath, target.name)
        align_res = ','.join(input_res_blocks)
        target.align_pdb = slice_out(inpath, prefix, target.native, align_res)
        if '-align_pdb' in target.extra_flags:
            target.extra_flags.pop('-align_pdb')
    elif '-align_pdb' in target.extra_flags:
        target.align_pdb = inpath+'/'+target.extra_flags['-align_pdb']
        assert( exists(align_pdb) )
    else:
        target.align_pdb = None

    # get sample loop res
    target.loop_res = {}

    if target.input_res == '-':
        if args.swa:
            print "WARNING: input_res[ name ] == '-' "
        continue

    ( workres , workchains  ) = parse_tag( target.working_res, alpha_sort=True )
    ( inputres , inputchains  ) = parse_tag( target.input_res, alpha_sort=True )

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
    ( workres , workchains  ) = parse_tag( target.working_res, alpha_sort=True )

    loopres_conventional = [ str(workchains[idx])+':'+str(workres[idx]) for idx in xrange( len( workres ) ) if (workres[idx] in loopres and workchains[idx] == loopchains[loopres.index(workres[idx])]) ]
    loopres_conventional = string.join( [ str(x) for x in loopres_conventional ] ,' ')
    target.loop_res[ 'conventional' ] = loopres_conventional

    if args.swa:
        loopres_swa = [ idx+1 for idx in xrange( len( workres ) ) if (workres[idx] in loopres and workchains[idx] == loopchains[loopres.index(workres[idx])]) ]
        loopres_swa = string.join( [ str(x) for x in loopres_swa ] ,' ')
        target.loop_res[ 'swa' ]  = loopres_swa


    # get VDW_rep_screen_info, it will only be used if -VDW_rep_screen_info flag is set in extra_flags_benchmark
    if '-VDW_rep_screen_info' in extra_flags_benchmark:

        periph_res_radius = 50.0
        if ( 'rrna' in target.name ) or ( 'rRNA' in target.name ):
            periph_res_radius = 100.0

        prefix = '%s/%s_%d_ANGSTROM_GRID_' % (inpath, target.name, periph_res_radius)
        target.VDW_rep_screen_pdb = prefix + target.native
        target.VDW_rep_screen_info = basename(target.VDW_rep_screen_pdb)

        if not exists( target.VDW_rep_screen_pdb ):

            loopres_list=string.split( target.loop_res[ 'conventional' ], ' ' )
            periph_res_tag = get_surrounding_res_tag( inpath+target.native, sample_res_list=loopres_list, radius=periph_res_radius, verbose=args.verbose )
            assert( len( periph_res_tag ) )
            slice_out( inpath, prefix, target.native, periph_res_tag )

            if args.verbose:
                print 'loopres_list for '+target.name+' = '+string.join(loopres_list)
                print 'periph_res for '+target.name+' = '+periph_res_tag



# write qsubMINIs, READMEs and SUBMITs
submit_files = init_submit_files()

for target in targets:
    
    dirname = target.name
    if not exists( dirname ): os.system( 'mkdir '+dirname )

    # move all required files to the correct directory
    start_files = target.helix_files + target.input_pdbs
    infiles = start_files + [ target.fasta, target.working_native ]
    if target.VDW_rep_screen_pdb:
        infiles.append(target.VDW_rep_screen_pdb)
    if target.align_pdb:
        if '-align_pdb' in extra_flags_benchmark:
            extra_flags_benchmark['-align_pdb'] = basename(target.align_pdb)
        infiles.append(target.align_pdb)
    if '-input_pdb' in target.extra_flags:
        input_pdb = inpath+'/'+target.extra_flags['-input_pdb']
        assert( exists(input_pdb) )
        infiles.append(input_pdb)
    os.system( 'cp %s %s/ ' % (' '.join(infiles), dirname) )

    # SETUP for StepWise Assembly
    if args.swa:

        fid = open( '%s/README_SWA' % dirname, 'w' )
        fid.write( os.environ['SWA_DAGMAN_TOOLS']+'/SWA_DAG/setup_SWA_RNA_dag_job_files.py' )
        if len( start_files ) > 0 :
            fid.write( ' -s' )
            for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
        if len( target.native ) > 0:
            fid.write( ' -native_pdb %s' % basename( target.working_native ) )
        fid.write( ' -fasta %s' %  basename( target.fasta ) )
        fid.write( ' -sample_res %s' % target.loop_res[ 'swa' ] )

        # case-specific extra flags
        for key, value in target.extra_flags.iteritems():
            flag = ' '.join([key, value]).strip()
            fid.write(' %s' % flag)

        # extra flags for whole benchmark
        for key, value in extra_flags_benchmark.iteritems():
            if '-analytic_etable_evaluation' in key:
                continue ### SWM Specific
            if '-motif_mode' in key:
                continue ### SWM Specific
            if '-score:weights' in key:
                key = '-force_field_file'
                weights_file = value
                if not exists( weights_file ):
                    weights_file = os.environ['ROSETTA_DB_WEIGHTS'] + weights_file
                assert( exists(weights_file) )
                os.system( 'cp %s %s' % (weights_file, target.name) )
            if '-score:rna_torsion_potential' in key:
                key = '-rna_torsion_potential_folder'
            if '-VDW_rep_screen_info' in key and 'True' in value:
                if target.VDW_rep_screen_info is None:
                    continue
                value = target.VDW_rep_screen_info
                fid.write(' -apply_VDW_rep_delete_matching_res False')
            flag = ' '.join([key, value]).strip()
            fid.write(' %s' % flag)
            
        fid.close()

        print '\nSetting up submission files for: ', target.name
        CWD = os.getcwd()
        fid_submit = open( dirname+'/SUBMIT_SWA', 'w' )
        fid_submit.write( os.environ['SWA_DAGMAN_TOOLS']+'/dagman/submit_DAG_job.py' )
        fid_submit.write( ' -master_wall_time %d' % args.nhours )
        fid_submit.write( ' -master_memory_reserve 2048' )
        fid_submit.write( ' -num_slave_nodes %d' % njobs )
        fid_submit.write( ' -dagman_file rna_build.dag' )
        fid_submit.close()

        for submit_file in submit_files:
            with open(submit_file,'a') as fid_submit:
                fid_submit.write( 'cd %s; source ./README_SWA && source ./SUBMIT_SWA; cd %s\n' % ( dirname, CWD ) )

    # SETUP for StepWise Monte Carlo
    else:

        fid = open( '%s/README_SWM' % target.name, 'w' )
        fid.write( os.environ['ROSETTA_BIN'] + 'stepwise @flags -out:file:silent swm_rebuild.out\n' )
        fid.close()

        fid = open( '%s/flags' % target.name, 'w' )
        if len( start_files ) > 0 :
            fid.write( '-s' )
            for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
            fid.write( '\n' )
        if len( target.native ) > 0:
            fid.write( '-native %s\n' % basename( target.working_native ) )
        if motif_mode_off:
            if len( target.terminal_res ) > 0:
                fid.write( '-terminal_res %s  \n' % make_tag_with_conventional_numbering( target.terminal_res, target.resnums, target.chains ) )
            if len( target.extra_min_res ) > 0 and not args.extra_min_res_off: ### Turn extra_min_res off for SWM when comparing to SWA
                fid.write( '-extra_min_res %s \n' % make_tag_with_conventional_numbering( target.extra_min_res, target.resnums, target.chains ) )
        #if ( len( input_pdbs[ name ] ) == 0 ):
        #    fid.write( '-superimpose_over_all\n' ) # RMSD over everything -- better test since helices are usually native
        fid.write( '-fasta %s\n' % basename( target.fasta) )
        if '-move' not in extra_flags_benchmark:
            if '-cycles' not in extra_flags_benchmark:
                fid.write( '-cycles 200\n' )
            if '-nstruct' not in extra_flags_benchmark:
                fid.write( '-nstruct 20\n' )
        #fid.write( '-intermolecular_frequency 0.0\n' )
        if not args.save_times_off:
            fid.write( '-save_times\n' )

        # case-specific extra flags
        for key, value in target.extra_flags.iteritems():
            flag = ' '.join([key, value]).strip()
            fid.write('%s\n' % flag)

        # extra flags for whole benchmark
        for key, value in extra_flags_benchmark.iteritems():
            if '-single_stranded_loop_mode' in key:
                continue ### SWA Specific
            if '-score:weights' in key:
                weights_file = value
                if not exists(weights_file):
                    weights_file = os.environ['ROSETTA_DB_WEIGHTS'] + weights_file
                assert( exists(weights_file) )
                os.system( 'cp %s %s' % (weights_file, target.name) )
            if '-VDW_rep_screen_info' in key and 'true' in value:
                if target.VDW_rep_screen_info is None:
                    continue
                value = basename(target.VDW_rep_screen_info)
            flag = ' '.join([key, value]).strip()
            fid.write('%s\n' % flag)
                
        fid.close()

        print '\nSetting up submission files for: ', target.name
        CWD = os.getcwd()
        os.chdir( target.name )

        rosetta_submit_cmd = 'rosetta_submit.py README_SWM SWM %d %d' % (njobs, args.nhours )
        if args.save_logs:
            rosetta_submit_cmd += ' -save_logs'
        os.system( rosetta_submit_cmd )

        os.chdir( CWD )

        for submit_file in submit_files:
            with open(submit_file,'a') as fid_submit:
                fid_submit.write( 'cd %s; source %s; cd %s\n' % ( target.name, submit_file,  CWD ) )


