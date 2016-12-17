#!/usr/bin/python

################################################################################
### import modules
################################################################################
import os
import operator
import argparse
import numpy as np
import subprocess as sp
import multiprocessing as mp
from os.path import exists, basename, dirname, isdir, abspath, expandvars
from glob import glob
from titles import get_title
from get_sequence import get_sequences


################################################################################
### GLOBAL DEFINITIONS
################################################################################
benchmark_dir = abspath(__file__).split('benchmark/')[0] + 'benchmark/'
column_labels = [
    " "*10,
    "Motif Properties",
    "Best of Five Lowest Energy Cluster Centers",
    "Lowest RMSD Model",
    "Lowest Energy Sampled" 
]
subcolumn_labels = [
    # " "*10 
    "Motif Name",    
    # Motif Properties
    "Length",
    "PDB",
    # Best of Five Lowest Energy Cluster Centers
    "Cluster Rank",
    "RMSD",
    "N_WC",
    "N_NWC",
    "fWC",
    "fNWC",
    "Rosetta Energy (RU)",
    # Lowest RMSD Model
    "RMSD", 
    "N_WC",
    "N_NWC",
    "fWC",
    "fNWC",
    # Lowest Energy Sampled
    "N_WC",
    "N_NWC",
    "fWC",
    "fNWC",
    "Rosetta Energy (RU)",
    "E-Gap to Opt. Exp. (RU)"
]


################################################################################
### CLASSES
################################################################################
class Command(object):

    def __init__(self, command, args=None):
        self.command = self._join(command, args)
        self.silent = False
        self._success = False
        self._out = None
        self._err = None
        self._stdout = sp.PIPE
        self._stderr = sp.PIPE
        self._outfile = None
        self._errfile = None

    def _join(self, x, y, delim=' '):
        if isinstance(x, list):
            x = delim.join(x)
        if isinstance(y, list):
            y = delim.join(y)
        if x is None and y is None:
            return None
        if x is None or y is None:
            return x if y is None else y
        return delim.join([str(x),str(y)]) 

    def _submit(self):
        pipe = sp.Popen( self.command, shell=True, 
                         stdout=self._stderr, 
                         stderr=self._stderr ) 
        self._out, self._err = pipe.communicate()
        self._check_error()
        self._save_logs()
        return

    def _check_error(self):
        self._success = True
        if self._err and len(self._err):
            if not self.silent:
                print "\nFAILED:", self.command
            self._success = False
        return self._success

    def _save_logs(self):
        if self._outfile:
            with open( self._outfile, 'w' ) as f:
                f.write( self._join(self.command, [self._out,self._err], delim='\n') )
        if self._errfile:
            with open( self._errfile, 'w' ) as f:
                if self._success is False:
                    f.write( self._join(self.command, self._err, delim='\n') )
        return

    def add_argument(self, argument, value=None):
        argument = self._join(argument, value)
        self.command = self._join(self.command, argument)
        return

    def save_logs(self, filename=None, out=True, err=True):
        if filename is None:
            filename = 'LOG_' + basename(self.command.split(' ')[0])
        filename = filename.split('.')[0]
        self._outfile = filename + '.out' if out is True else None
        self._errfile = filename + '.err' if err is True else None
        return

    def submit(self):
        self._submit()
        return self._success
    
    def output(self):
        self._submit()
        return self._out if self._success else None

################################################################################
class TableRow(object):

    def __init__(self):
        self._columns = []

    def add_columns(self, cols):
        if not isinstance(cols, list):
            cols = [ cols ]
        for col in cols:
            self._columns.append( col )
        return

    def columns(self):
        return self._columns

################################################################################
class Table(object):

    def __init__(self, filename):
        self._rows = []
        self._data_rows = []
        self.filename = filename
        self.delimiter = '\t'
        self.newline = '\n'

    def _format_string(self, value, width=None):
        if width is None:
            width = len(str(value))
        if isinstance(value, int):
            return "%-*d" % (width, value)
        if isinstance(value, float):
            return "%-*.2f" % (width, value)
        value = "--" if value is None else str(value)
        return "%-*s" % (width, str(value))
        
    def _row_to_string(self, row):
        return self.delimiter.join(map(self._format_string, row)) 
    
    def _table_to_string(self):
        return self.newline.join(map(self._row_to_string, self._rows))

    def add_row(self, row):
        if not isinstance(row, list):
            row = [ row ]
        self._rows.append( row )
        return

    def add_data_row(self, row):
        if not isinstance(row, list):
            row = [ row ]
        self._data_rows.append( row )
        self._rows.append( row )
        return

    def column_averages(self):
        col_accums = []
        col_counts = []
        for data_idx, data in enumerate(self._data_rows):
            idx = 0
            for col in data:
                if col is None or isinstance(col, str):
                    increment = 0.0
                else:
                    col = float( col )
                    increment = 1.0
                if data_idx == 0:
                    col_accums.append( col )
                    col_counts.append( increment )
                else:
                    col_accums[idx] += col
                    col_counts[idx] += increment
                idx += 1
        col_aves = []
        for accum, count in zip(col_accums,col_counts):
            if count == 0.0:
                col_aves.append(None)
                continue
            col_aves.append(accum/count)
        return col_aves

    def save(self):
        with open( self.filename, 'w' ) as f:
            f.write( self._table_to_string() )
        return 

    def merge_tables(self, filenames):
        if not isinstance(filenames, list):
            filenames = [ filenames ]
        subtable_rows = []
        for table_idx, filename in enumerate(filenames):
            with open( filename, 'r' ) as fin:
                 for idx, line in enumerate(fin.readlines()):
                    cols = line.strip().split(self.delimiter)
                    if table_idx < 1:
                        subtable_rows.append(cols)
                    else:
                        startcol = 1 if idx == 0 else 3
                        if cols[0] != subtable_rows[idx][0]:
                            continue
                        subtable_rows[idx] += cols[startcol:]
        self.add_row( [basename(f).split('.')[0] for f in filenames] )
        for row in subtable_rows:
            self.add_row( row )
        return


################################################################################
### HELPER FUNCTIONS
################################################################################
def find( file, start='./' ):
    start += '/'
    if exists(start+file):
        return start+file
    for dir in filter(isdir, glob(start+'*')):
        return find(file, start=dir)
    return None

################################################################################
def get_rosetta_exe( exe, tools=False ):
    rosetta = expandvars('$ROSETTA')
    if '$' in rosetta:
        rosetta = '~/src/rosetta/'
    if tools:
        rosetta += '/tools/'
    else:
        rosetta += '/main/source/bin/'
    return find( exe, start=rosetta )

################################################################################
def get_target_names():
    target_names = []
    info_files = glob(benchmark_dir+'/input_files/*.txt')
    for file in info_files:
        with open( file, 'r' ) as f:
            for line in f:
                if 'Name:' not in line:
                    continue
                name = line.split()[-1].strip()
                if name in target_names:
                    continue
                target_names.append( name )
    return target_names 

################################################################################
def get_score_data( filename, colnames=['score'], sort=None, filters=None, tags=None, keep=None ):
    """
    Gets a bunch of data out of a silent file. If the data is not present in the silent file,
    opens the base pair analysis version fo that silent file. If that doesn't exist, make it.
    """

    # UGH: We are going to have to do two separate read-throughs then join arrays
    # indexed by the same key into one.
    # Wait, actually... let's assume the NEW ENERGIES are going to be good. May
    # limit legal SHAs.
    # AMW: I have shunted this off onto build_full_model. I am a hero.
    bps_silent_file = filename #create_bps_silent_file( filename )
    
    if not isinstance(colnames, list):
        colnames = [ colnames ]
    data = []
    colidx = None
    with open( bps_silent_file, 'r' ) as f:
        for line in f:
            if not "SCORE:" in line:
                continue
            cols = filter(None,[c.strip() for c in line.split()])
            if not len(cols):
                continue
            if "description" in line:
                colidx = map(cols.index, filter(cols.count, colnames))
                continue
            if tags and not any( t in line for t in tags ):
                continue
            data.append([])
            for idx in colidx:
                col = cols[idx]
                try:
                    data[-1].append( float(col) )
                except:
                    data[-1].append( col )
            data[-1] = tuple(data[-1])
    
    if filters is not None:
        if not isinstance(filters, list):
            filters = [ filters ] 
        for idx, value in enumerate(filters):
            if value is None:
                continue
            data = [d for d in data if d[idx] <= value]
    if not len(data):
        return None
    if sort is not None:
        sort = colnames.index(sort) if isinstance(sort, str) else sort-1
        data = sorted(data, key=operator.itemgetter(sort))
    if len(colidx) == 1:
        data = [d[0] for d in data] 
    if keep is not None:
        keep = min(keep, len(data))
        data = data[0] if keep == 1 else data[:keep]
    return data

################################################################################
def get_rmsd_type( silent_file ):
    with open( silent_file, 'r' ) as f:
        for line in f:
            if not 'SCORE:' in line:
                continue
            if 'rms_fill' in line:
                return 'rms_fill'
            elif 'NAT_rmsd' in line:
                return 'NAT_rmsd'
            else:
                break
    return 'rms'

################################################################################
def get_working_target():
    return basename(os.getcwd())

################################################################################
def get_silent_file( filename=['region_FINAL.out','swm_rebuild.out'], dir=None ):
    if not isinstance(filename, list):
        filename = [ filename ]
    if dir is not None:
        filename = ['/'.join([dir,file]) for file in filename]
    silent_files = filter(exists, filename)
    return silent_files[0] if len(silent_files) else None

################################################################################
def get_native_pdb():
    target = get_working_target()
    try:
        native_pdb = glob( target+'_????_RNA.pdb' )[0]
    except:
        try:
            native_pdb = glob( target+'_NATIVE_????_RNA.pdb' )[0]
        except:
            try:
                native_pdb = glob( target+'_NATIVE_*.pdb' )[0]
            except:
                print target
    return native_pdb

################################################################################
def get_start_pdb_list():
    target = get_working_target()
    try:
        start_pdbs = glob( target+'_START*_????_RNA.pdb' )
    except:
        start_pdbs = glob( target+'_START*.pdb' )
    return start_pdbs

################################################################################
def get_motif_length():
    length = len(''.join(get_sequences(get_native_pdb())[0]))
    for start_pdb in get_start_pdb_list():
        length -= len(''.join(get_sequences(start_pdb)[0]))
    print length
    return length

################################################################################
def get_pdb_id():
    print get_native_pdb
    return get_native_pdb().split('_')[-2].upper()

################################################################################
def get_opt_exp_score( inpaths ):
    target = get_working_target()
    opt_exp_score = None
    for inpath in inpaths:
        inpath = find( '/'.join([inpath,target]), start='../../' )
        if inpath is None:
            continue
        silent_file = get_silent_file( dir=inpath )
        if silent_file is None:
            continue
        score_types = ['score', get_rmsd_type(silent_file)]
        cutoffs = [None, 1.5]
        data = get_score_data( silent_file, colnames=score_types, sort='score', filters=cutoffs, keep=1 )
        if data is None:
            continue 
        score = data[score_types.index('score')]
        if opt_exp_score is not None and score >= opt_exp_score:
            continue
        opt_exp_score = score
    return opt_exp_score

################################################################################
def get_flag( flag ):
    if exists( 'flags' ):
        with open( 'flags', 'r' ) as f:
            for line in f:
                if flag not in line: 
                    continue
                return line.strip()
    if exists( 'README_SWA' ):
        with open( 'README_SWA', 'r' ) as f:
            for line in f.read().strip().split(' -'):
                line = '-'+line
                if flag not in line:
                    continue
                return line.strip()
    return None     

################################################################################
def virtualize_missing_residues( silent_file ):
    silent_file_out = silent_file.replace(".out","_full_model.out")
    build_full_model_exe = get_rosetta_exe( "build_full_model" )
    weights = None #get_flag( "-score:weights" ).split(' ')[-1]
    torsion_potential = None #get_flag( "-score:rna_torsion_potential" ).split(' ')[-1]
    command = Command( build_full_model_exe )
    command.add_argument( "-in:file:silent", value=silent_file )
    command.add_argument( "-out:file:silent", value=silent_file_out )
    command.add_argument( "-in:file:native", value=get_native_pdb() ) # for native base pairs
    command.add_argument( "-out:overwrite", value="true" )
    if weights is not None:
        command.add_argument( "-score:weights", value=weights )
    if torsion_potential is not None:
        command.add_argument( "-score:rna_torsion_potential", value=torsion_potential )
    command.add_argument( "-virtualize_built", value="true" )
    command.save_logs()
    success = command.submit()
    return silent_file_out if success is True else None

################################################################################
def get_full_model_parameter(silent_file, parameter):
    with open( silent_file, 'r' ) as f:
        for line in f:
            if "FULL_MODEL_PARAMETERS" not in line:
                continue
            if parameter not in line:
                continue
            cols = line.strip().split()
            return cols[cols.index(parameter)+1]         
    return None

################################################################################
def make_res_list(tag):
    new_list = []
    subtags = tag.replace(',',' ').split(' ')
    for subtag in subtags:
        if '-' in subtag:
            subtag = subtag.split('-')
            for x in xrange(int(subtag[0]),int(subtag[1])+1):
                new_list.append(str(x))
        else:
            res_list.append(subres)
    return new_list

################################################################################
def make_res_tag(res, exclude=None, delim=' ', dashes=True):
    exclude_list = make_res_list(exclude) if exclude else []
    res_list = [x for x in make_res_list(res) if x not in exclude_list]
    if dashes:
        for idx, res in enumerate(res_list):
            while idx+2 != len(res_list) and int(res_list[idx+1]) + 1 == int(res_list[idx+2]):
                res_list.pop(idx+1)
            res_list[idx] = res_list[idx] + '-' + res_list.pop(idx+1)
    res_tag = delim.join(res_list)
    return res_tag
    
################################################################################
def create_common_args_file( silent_file ):
    common_args_file = 'SWA_cluster_common_args_'+silent_file
    if exists( common_args_file ):
        Command( "rm -f ", args=common_args_file ).submit()
    working_res = get_full_model_parameter(silent_file, 'WORKING')
    sample_res = get_full_model_parameter(silent_file, 'CALC_RMS')
    if sample_res is None:
        sample_res = get_full_model_parameter(silent_file, 'SAMPLE')
    fixed_res = make_res_tag(working_res, exclude=sample_res)
    common_args = []
    common_args.append( '-in:file:silent_struct_type binary_rna' )
    common_args.append( get_flag('-score:weights') )
    common_args.append( get_flag('-score:rna_torsion_potential') )
    common_args.append( get_flag('-fasta') )
    if get_flag('-VDW_rep_screen_info'):
        common_args.append( get_flag('-VDW_rep_screen_info') )
        common_args.append( '-VDW_rep_delete_matching_res false' )
    common_args.append( '-jump_point_pairs ' + working_res )
    common_args.append( '-alignment_res ' + working_res ) 
    common_args.append( '-fixed_res ' + fixed_res )
    common_args.append( '-input_res ' + fixed_res )
    common_args.append( '-input_res2 ' + sample_res )
    common_args.append( '-global_sample_res_list ' + sample_res ) 
    common_args.append( '-sample_res ' + sample_res ) 
    common_args.append( '-rmsd_res ' + sample_res ) 
    with open( common_args_file, 'w' ) as fid:
        fid.write(' ')
        fid.write(' '.join(filter(None, common_args)))
        fid.write('\n')
    if not exists( common_args_file ):
        common_args_file = None
    return common_args_file

def create_bps_silent_file( silent_file ):
    # todo: don't run on FARFAR because it already has it.
    analysis_exe = get_rosetta_exe( 'analyze_base_pairing' )
    analysis_silent_file = silent_file.replace( '.out', '_base_pairing_added.out' )
    # THIS IS EXPENSIVE. Don't remake, especially not in these exploratory
    # investgations...
    #if exists( analysis_silent_file ):
    #    Command( "rm -f ", args=analysis_silent_file ).submit()
    command = Command( analysis_exe )
    command.add_argument( "-in:file:silent", value=silent_file ) 
    command.add_argument( "-out:file:silent", value=analysis_silent_file ) 
    command.save_logs()
    success = command.submit()
    return analysis_silent_file if success is True else None

################################################################################
def create_cluster_silent_file( silent_file ):
    common_args_file = create_common_args_file( silent_file )
    if 'swm' in silent_file:
        silent_file_virt = virtualize_missing_residues( silent_file )
        if not silent_file_virt or not exists( silent_file_virt ):
            return silent_file
        silent_file = silent_file_virt
    cluster_rmsd = 2.0 
    suite_cluster_rmsd = 2.5 
    rename_tags = False
    no_graphic = False
    ignore_unmatched_virtual_res = False
    native_pdb = get_native_pdb()
    cluster_exe = "SWA_RNA_python/SWA_dagman_python/misc/SWA_cluster.py"
    cluster_exe = get_rosetta_exe( cluster_exe, tools=True )
    top_energy_clusters_folder = "TOP_ENERGY_CLUSTERS/"
    cluster_silent_file = "%s/top_energy_clusters.out" % top_energy_clusters_folder
    if exists( cluster_silent_file ):
        Command( "rm -f ", args=cluster_silent_file ).submit()
    elif not exists( top_energy_clusters_folder ):
        Command( "mkdir -p", args=top_energy_clusters_folder ).submit()
    command = Command( cluster_exe )
    command.add_argument( "-num_pose_kept", value=100 )
    command.add_argument( "-distinguish_pucker", value="false" )
    command.add_argument( "-extract_pdb ", value="False" ) 
    command.add_argument( "-cluster_rmsd", value=cluster_rmsd )
    command.add_argument( "-suite_cluster_rmsd", value=suite_cluster_rmsd )
    command.add_argument( "-silent_file", value=silent_file )
    command.add_argument( "-native_pdb", value=native_pdb )
    command.add_argument( "-output_filename", value=cluster_silent_file )    
    command.add_argument( "-full_length_loop_rmsd_clustering", value="True" )
    if not rename_tags:
        command.add_argument( "-clusterer_rename_tags", value="false" )
        command.add_argument( "-add_lead_zero_to_tag", value="false" )
    if not no_graphic:
        command.add_argument( "-no_graphic", value="False" )
    if ignore_unmatched_virtual_res:
        command.add_argument( "-ignore_unmatched_virtual_res", value="True" )
    if common_args_file:
        command.add_argument( "-common_args", value=common_args_file )
    command.save_logs()
    success = command.submit()
    return cluster_silent_file if success is True else None

################################################################################
def get_lowest_energy_cluster_centers( nclusters=5 ):
    silent_file = get_silent_file()    
    cluster_silent_file = create_cluster_silent_file( silent_file )
    if cluster_silent_file is None or not exists( cluster_silent_file ):
        print "\n[WARNING] cluster_silent_file not found for target: %s" % get_working_target()
        cluster_silent_file = silent_file
    cluster_center_list = []
    # AMW: plan is for get_score_data to know where to look.
    score_types = ['score', get_rmsd_type(cluster_silent_file), 'N_WC', 'N_NWC', 'f_natWC', 'f_natNWC' ]
    data = None
    if 'swm' in silent_file:
        # PROBLEM: we need to build/virtualize missing residues in SWM silent files before 
        # clustering with SWA_cluster.py, but the scores before and after build_full_model 
        # are not matching up (yet!) ... 
        # SOLUTION: use 'build_full_model -virtualize_built' output for clustering, but get
        # original scores of poses, using the tags found in clusterer output
        tags = get_score_data( cluster_silent_file, colnames='description' )
        data = get_score_data( silent_file, colnames=score_types, sort='score', tags=tags, keep=nclusters )
    else:
        data = get_score_data( cluster_silent_file, colnames=score_types, sort='score', keep=nclusters )
    if data is None:
        return cluster_center_list
    for idx, (energy, rmsd, nwc, nnwc, fwc, fnwc ) in enumerate(data, start=1):
        if idx > nclusters:
            break
        # Following order in subcolumns
        cluster_center_list.append( [idx, rmsd, nwc, nnwc, fwc, fnwc, energy] )
    return cluster_center_list


################################################################################
### MAIN FUNCTIONS
################################################################################
def get_target_properties():
    ''' 
        column:    | Motif Properties |
        subcolumn: | Length | PDB     |

    '''
    length = get_motif_length()
    pdb_id = get_pdb_id()
    return [ length, pdb_id ]

################################################################################
def get_best_of_lowest_energy_cluster_centers():
    '''
        column:    | Best of Five Lowest Energy Cluster Centers |
        subcolumn: | Cluster Rank | RMSD | N_WC | N_NWC | fWC | fNWC | Rosetta Energy (RU)  |

    '''
    cluster_centers = get_lowest_energy_cluster_centers()
    for idx, cluster_center in enumerate(cluster_centers):
        rmsd = cluster_center[1]
        if idx and rmsd >= best_rmsd:
            continue
        cluster_centers[0] = cluster_center
        best_rmsd = rmsd
    return cluster_centers[0]

################################################################################
def get_lowest_rmsd_model():
    '''
        column:    | Lowest RMSD Model |
        subcolumn: | RMSD              | N_WC | N_NWC | fWC | fNWC |

    '''
    silent_file = get_silent_file()
    rmsd_type = get_rmsd_type( silent_file )
    scoretypes = [ rmsd_type, "N_WC", "N_NWC", "f_natWC", "f_natNWC" ]
    rmsd = get_score_data( silent_file, colnames=scoretypes, sort=rmsd_type, keep=1 )
    return [ rmsd ]

################################################################################
def get_lowest_energy_sampled( opt_exp_inpaths ):
    '''
        column:    | Lowest Energy Sampled                         |
        subcolumn: | N_WC | N_NWC | fWC | fNWC | Rosetta Energy (RU) | E-Gap to Opt. Exp. (RU) |

    '''
    silent_file = get_silent_file()
    scoretypes = [ "score", "N_WC", "N_NWC", "f_natWC", "f_natNWC" ]
    energy, nwc, nnwc, fwc, fnwc = get_score_data( silent_file, colnames=scoretypes, sort='score', keep=1 )
    if energy is None:
        return [ None, None, None, None, None, None ]
    opt_exp_energy = get_opt_exp_score( opt_exp_inpaths )
    if opt_exp_energy is None:
        return [ energy, nwc, nnwc, fwc, fnwc, None ]
    energy_gap = energy - opt_exp_energy
    return [ energy, nwc, nnwc, fwc, fnwc, energy_gap ]


