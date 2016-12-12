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
parser.add_argument('-j','--njobs', default='10', type=int, help='Number of cores for each job.')
parser.add_argument('--swa', action='store_true', help='Additional flag for setting up SWA runs.')
parser.add_argument('--farna', action='store_true', help='Additional flag for setting up FARNA runs (no minimize).')
parser.add_argument('--farfar', action='store_true', help='Additional flag for setting up FARFAR runs (FARNA+minimize).')
parser.add_argument('--extra_min_res_off', action='store_true', help='Additional flag for turning extra_min_res off.')
parser.add_argument('--save_times_off', action='store_true', help='Additional flag for turning save_times flag off.')
parser.add_argument('-motif_mode_off', help="temporary hack for turning off hardcoded '-motif_mode' flag", action="store_true")
#parser.add_argument('--path_to_rosetta', default='', help='Path to working copy of rosetta.')
parser.add_argument('--rosetta', default='', help='Path to working copy of rosetta.')
parser.add_argument('--no_align_pdb', action='store_true', help='Do not read in align_pdb; forces alignment to native, useful for native_screen runs.')
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
parser.add_argument('--save_logs', help="save .out and .err logs for each job.", action="store_true")
parser.add_argument('-queue', default='', help='Queue for cluster submission.')
parser.add_argument('--block_stack_off', action='store_true', help='Additional flag for turning off -block_stack_above_res & -block_stack_below_res.')
parser.add_argument('-stepwise_lores',action='store_true', help='used to setup stepwise_lores mode (FARNA optimization)')
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
if not exists( info_file ): info_file = dirname(argv[ 0 ]) + "/../../../input_files/" + info_file
print info_file
assert( exists( info_file ) )


# define and check paths
inpath = info_file.replace('.txt', '' ) + '/'
assert( exists( inpath ) )

def onelettersequence( sequence ):
    # replace nonstandard names like Z[Mg] with Z
    sequence_clean = []
    in_non_standard = False
    for c in sequence:
        if c == '[':
            in_non_standard = True
            continue
        if c == ']':
            in_non_standard = False
            continue
        if not in_non_standard: sequence_clean += c
    return sequence_clean

# read info_file
info_fid = file_handlers.TargetDefinitionsFile()
info_fid.load(open(info_file))
assert( info_fid.validate() )
target_definitions = info_fid.target_definitions
if len(args.user_input_runs):
    target_definitions = [td for td in target_definitions if td.name not in args.user_input_runs]
targets = [info_handlers.Target(td) for td in target_definitions]

# If we weren't given extra_flags_benchmark.txt on the command line
# use the "new format file" 'Benchmark_flags:' cue instead.
if extra_flags_benchmark == {}:
    extra_flags_benchmark = info_fid.extra_flags_benchmark
if input_res_benchmark is None:
    input_res_benchmark = info_fid.input_res_benchmark
if extra_min_res_benchmark is None:
    extra_min_res_benchmark = info_fid.extra_min_res_benchmark


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
    # We need to sort the sequences first. That's absurd, of course, but otherwise we'd need a 
    # vastly more clever sequence determination algorithm. That's a TODO.
    assert( sorted(string.join(sequences,'')) == sorted(string.join(get_sequences( target.working_native )[0],'')) )

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
    # This is going to be a huge pain to handle so just ignore
    #assert( chainbreak_pos == chainbreak_pos_secstruct )
    stems = get_all_stems( secstruct_joined, chainbreak_pos, sequence_joined  )

    # instead of sequence_joined, we have to look at fasta_entities to get the real story,
    # NCNTs and all
    fasta_entities = []
    i = 0
    while i < len(sequence_joined):
        if i == len(sequence_joined) - 1 or sequence_joined[i+1] != '[':
            # simple case
            fasta_entities.append(sequence_joined[i])
        else:
            entity = sequence_joined[i]
            i += 1
            while sequence_joined[i] != ']':
                entity += sequence_joined[i]
                i += 1
                entity += ']'
            fasta_entities.append(entity)
        i += 1

    for i in range( len( stems ) ):
        helix_file =  '%s/%s_HELIX%d.pdb' % (inpath,target.name,(i+1))

        stem = stems[i]
        print stems[i]
        helix_seq = ''; helix_resnum = [];
        for bp in stem:
            helix_seq    += fasta_entities[ bp[0] - 1 ] #sequence_joined[ bp[0] - 1 ]
            helix_resnum.append( bp[0] )
        helix_seq += ' '
        for bp in stem[::-1]:
            helix_seq    += fasta_entities[ bp[1] - 1 ]
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
    # AMW: recall that the length of the sequence is probably a lie.
    #L = len( sequence_joined )
    L = len( fasta_entities )
    target.terminal_res = []
    target.extra_min_res = []
    target.block_stack_above_res = []
    target.block_stack_below_res = []
    for m in range( 1, L+1 ):
        if ( m not in input_resnum_fullmodel ): continue
        prev_moving = ( m - 1 not in input_resnum_fullmodel ) and ( m != 1 )
        next_moving = ( m + 1 not in input_resnum_fullmodel ) and ( m != L )
        right_before_chainbreak = ( m == L or m in chainbreak_pos )
        right_after_chainbreak  = ( m == 1 or m - 1 in chainbreak_pos )
        if ( ( right_after_chainbreak and not next_moving ) or \
            ( right_before_chainbreak and not prev_moving ) ):
            target.terminal_res.append( m )
            if right_after_chainbreak:
                target.block_stack_below_res.append( m )
            if right_before_chainbreak:
                target.block_stack_above_res.append( m )
        if ( ( prev_moving and not next_moving and not right_before_chainbreak ) or \
            ( next_moving and not prev_moving and not right_after_chainbreak ) ):
            target.extra_min_res.append( m )

    if not '-motif_mode\n' in extra_flags_benchmark.keys() and not motif_mode_off:
        extra_flags_benchmark[ '-motif_mode' ] = ''# \n' )

    # needed for stepwise_lores to work with base pair steps that include flanking helices:
    if args.stepwise_lores:
        jump_res[ name ] = []
        cutpoint_closed[ name ] = []
        jump_bps = []
        stems = get_all_stems( secstruct_gen_joined )
        if len( stems ) == 0: stems = get_all_stems( secstruct_joined )
        for i in range( len( stems ) ):
            stem = stems[i]
            for bp in stem:
                jump_bps.append( bp )
                jump_res[ name ].extend( [ bp[ 0 ], bp[ 1 ] ] )
                if ( bp != stem[ -1] ): cutpoint_closed[ name ].append( bp[ 0 ] )
        # need to be really explicit about cutpoints_closed in stepwise right now... there is an edge case (srl_fixed)
        # where jumps don't quite work.
        for i in range( len( input_resnums_by_block ) ):
            cuts  = []
            for m in range(1, len( input_resnums_by_block[i] ) ):
                m_full = get_fullmodel_number( (input_resnums_by_block[i][m-1],input_chains_by_block[i][m-1]),resnums[name],chains[name] )
                if ( ( input_resnums_by_block[ i ][ m ] != input_resnums_by_block[ i ][ m-1 ] + 1 ) or
                     ( input_chains_by_block[ i ][ m ]  != input_chains_by_block[ i ][ m-1 ] ) or
                     ( m_full in cutpoint_closed[ name ] ) ):
                    cuts.append( m_full )
            for jump_bp in jump_bps:
                if ( ( target.resnums[ jump_bp[0]-1 ], target.chains[ jump_bp[0]-1 ] ) in zip( input_resnums_by_block[i], input_chains_by_block[i] ) and \
                     ( target.resnums[ jump_bp[1]-1 ], target.chains[ jump_bp[1]-1 ] ) in zip( input_resnums_by_block[i], input_chains_by_block[i] ) ):
                    cut_exists_for_jump = False
                    for cut in cuts:
                        if ( cut >= jump_bp[0] and cut < jump_bp[1] ): cut_exists_for_jump = True
                    if not cut_exists_for_jump:
                        cutpoint_closed[ name ].append( jump_bp[ 0 ] )
                        cuts.append( jump_bp[ 0 ] )

    # for cases that require docking two chains & farna, need to recognize and set up chain_connections flag
    target.dock_partners = []
    if args.farna or args.farfar:
        strand_num = 1
        strand = {}
        for m in range( 1, L+1 ):
            strand[m] = strand_num
            if ( m in chainbreak_pos ): strand_num += 1
        # which strands are connected up?
        strand_connected = {}
        for m in range(1,strand_num+1):
            strand_connected[ m ] = {}
            for n in range(1,strand_num+1): strand_connected[ m ][ n ] = False
        for domain in input_resnum_fullmodel_by_block:
            for m in domain:
                for n in domain:
                    strand_connected[ strand[m] ][ strand[n] ] = True
        # are there any separated clusters?
        already_in_cluster = {}
        clusters = []
        for m in range(1,strand_num+1): already_in_cluster[ m ] = False
        for m in range(1,strand_num+1):
            if already_in_cluster[ m ]: continue
            strand_cluster = set( [ m ] )
            for n in range(1,strand_num+1):
                if strand_connected[ m ][ n ]:
                    strand_cluster.add( n )
                    already_in_cluster[ n ] = True
            clusters.append( strand_cluster )
        assert( len( clusters ) > 0 )
        if len( clusters ) > 1:
            assert( len( clusters ) == 2 )
            dock_partners[ name ] = [ [], [] ]
            for m in range( 1, L+1 ):
                if strand[ m ] in clusters[ 0 ]:
                    dock_partners[ name ][ 0 ].append( m )
                else:
                    assert( strand[ m ] in clusters[ 1 ])
                    dock_partners[ name ][ 1 ].append( m )

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
        # For some reason, this doesn't try to create the file for you. I defy this and have added the next 3 lines.
        #prefix = '%s/%s_ALIGN_' % (inpath, target.name)
        #align_res = ','.join(input_res_blocks)
        #target.align_pdb = slice_out(inpath, prefix, target.native, align_res)
        target.align_pdb = inpath+'/'+target.extra_flags['-align_pdb']
        assert( exists(target.align_pdb) )
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

    def add_block_stack_flags( args, extra_flags, block_stack_above_res, block_stack_below_res, fid ):
        # used in FARNA & SWM
        if not args.block_stack_off and extra_flags[ name ].find( '-block_stack_off') == -1:
            if len( block_stack_above_res[ name ] ) > 0:
                fid.write( '-block_stack_above_res %s  \n' % make_tag_with_conventional_numbering( block_stack_above_res[ name ], resnums[ name ], chains[ name ] ) )
            if len( block_stack_below_res[ name ] ) > 0:
                fid.write( '-block_stack_below_res %s  \n' % make_tag_with_conventional_numbering( block_stack_below_res[ name ], resnums[ name ], chains[ name ] ) )
        return

    def add_start_files_flag( fid, start_files ):
        if len( start_files ) > 0 :
            fid.write( '-s' )
            for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
            fid.write( '\n' )


    def copy_extra_files( tag, cols, name ):
        # This had assumed that we are passing a vector not a dict as cols. Huh.
        if tag in cols.keys():
            #filename = cols[ cols.index( tag )+1 ]
            filename = cols[ tag ]
            print filename
            assert( exists( inpath+'/'+filename ) )
            system( 'cp %s/%s %s' % (inpath, filename, name ) )

    def add_extra_flags_for_name( fid, extra_flags, name ):
        # case-specific extra flags
        if len( extra_flags ) == 0 or extra_flags == '-': return

        #fid.write( '%s\n' % extra_flags[name] )
        #cols = extra_flags.split( ' ' )

        if not args.no_align_pdb: copy_extra_files( '-align_pdb', extra_flags, name )
        copy_extra_files( '-extra_res_fa', extra_flags, name)

        # These are already all parsed.
        for flag in extra_flags.keys():#parse_flags( extra_flags ):
            flag = flag.replace('True','true').replace('False','false')
            if flag == '-block_stack_off\n': continue
            if ( args.no_align_pdb and flag.find( '-align_pdb' ) > -1 ): continue
            fid.write( flag )

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

    elif args.farna or args.farfar: # Fragment Assembly of RNA

        fid = open( '%s/README_FARFAR' % name, 'w' )
        fid.write( 'rna_denovo @flags -out:file:silent farna_rebuild.out\n' )
        fid.close()

        fid = open( '%s/flags' % name, 'w' )
        fid.write( '-fasta %s.fasta\n' % name )
        if len( native[ name ] ) > 0:
            fid.write( '-working_native %s\n' % basename( working_native[ name ] ) );
        add_start_files_flag( fid, start_files )
        fid.write( '-working_res %s\n' % working_res[ name ].replace( ',',' ') )
        add_block_stack_flags( args, extra_flags, block_stack_above_res, block_stack_below_res, fid )
        if len( extra_min_res[ name ] ) > 0 and not args.extra_min_res_off:
            fid.write( ' -extra_minimize_res %s\\\n' % make_tag_with_conventional_numbering( extra_min_res[ name ], target.resnums, target.chains ) )
        if args.farfar: fid.write( '-minimize_rna true\n' )
        else: fid.write( '-minimize_rna false\n' )
        fid.write( ' -no_minimize\\\n' )
        if not cycles_flag_found:   fid.write( ' -cycles 20000\\\n' )
        if not nstruct_flag_found:  fid.write( ' -nstruct 20\\\n' )
        if not args.save_times_off: fid.write( ' -save_times\\\n' )
        # case-specific extra flags
        if ( len( extra_flags[name] ) > 0 ) and ( extra_flags[ name ] != '-' ) :
            for flag in parse_flags_string( extra_flags[ name ] ):
                flag = flag.replace('true','True').replace('false','False')
                fid.write( ' %s\\\n' % flag[:-1] )
        # silly, currently required for FARNA, but hopefully not in future
        #fid.write( '-output_res_num %s\n' % make_tag_with_dashes( target.resnums, target.chains ) )
        if len( target.dock_partners ) > 0:
            fid.write( '-chain_connection SET1 %s SET2 %s\n' % ( make_tag_with_conventional_numbering( target.dock_partners[ 0 ], target.resnums, target.chains ), \
                make_tag_with_conventional_numbering( target.dock_partners[ 1 ], target.resnums, target.chains ) ) )
                                                        
        add_extra_flags_for_name(fid, target.extra_flags, target.name)

        # extra flags for whole benchmark
        for flag in extra_flags_benchmark:
            if ( '-motif_mode' in flag ): continue # not really checked in FARFAR; use block_stack + extra_min_res instead.
            if ( '#' in flag ): continue
            flag = flag.replace('True','true').replace('False','false')
            fid.write( flag[:-1]+'\n' )
        fid.close()

        print '\nSetting up submission files for: ', name
        CWD = getcwd()
        chdir( name )

        rosetta_submit_cmd = 'rosetta_submit.py README_FARFAR FARFAR %d %d' % (njobs, args.nhours )
        if args.save_logs: rosetta_submit_cmd += ' -save_logs'
        if len(args.queue) > 0: rosetta_submit_cmd += ' -queue %s' % args.queue
        system( rosetta_submit_cmd )

        chdir( CWD )

        fid_qsub.write( 'cd %s; source %s; cd %s\n' % ( name, qsub_file,  CWD ) )
        if len( qsub_file2 ) > 0: fid_qsub2.write( 'cd %s; source %s; cd %s\n' % ( name, qsub_file2,  CWD ) )

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
        # Copied from master.
        if not args.block_stack_off and '-block_stack_off' not in extra_flags_benchmark:
            if len( target.block_stack_above_res ) > 0:
                fid.write( '-block_stack_above_res %s  \n' % make_tag_with_conventional_numbering( target.block_stack_above_res, target.resnums, target.chains ) )
            if len( target.block_stack_below_res ) > 0:
                fid.write( '-block_stack_below_res %s  \n' % make_tag_with_conventional_numbering( target.block_stack_below_res, target.resnums, target.chains ) )
        #add_start_files_flag( fid, start_files )
        
        if args.stepwise_lores:
            if ( len( jump_res[ name ] ) > 0 ):
                fid.write( '-jump_res %s \n' % make_tag_with_conventional_numbering( target.jump_res, target.resnums, target.chains ) )
            if ( len( cutpoint_closed[ name ] ) > 0 ):
                fid.write( '-cutpoint_closed %s \n' % make_tag_with_conventional_numbering( target.cutpoint_closed, target.resnums, target.chains ) )
            fid.write( '-include_neighbor_base_stacks\n' ) # Need to match FARNA.

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
      
        add_extra_flags_for_name(fid, target.extra_flags, target.name)

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

