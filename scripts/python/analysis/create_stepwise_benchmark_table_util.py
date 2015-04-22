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
	"Rosetta Energy (RU)",
	# Lowest RMSD Model
	"RMSD", 
	# Lowest Energy Sampled
	"Rosetta Energy (RU)",
	"E-Gap to Opt. Exp. (RU)"
]


################################################################################
### CLASSES
################################################################################
class Command(object):

	def __init__(self, command, args=None):
		self.command = self._join(command, args)
		self.stdout = sp.PIPE
		self.stderr = sp.PIPE
		self.silent = False
		self._out = None
		self._err = None

	def _join(self, x, y):
		if isinstance(x, list):
			x = ' '.join(x)
		if isinstance(y, list):
			y = ' '.join(y)
		return x if y is None else ' '.join([str(x),str(y)]) 

	def _run(self):
		pipe = sp.Popen( self.command, shell=True, 
						 stdout=self.stderr, 
						 stderr=self.stderr ) 
		self._out, self._err = pipe.communicate()
		return

	def _check_error(self):
		if self._err and len(self._err):
			if not self.silent:
				print '\n'.join([self.command,self._out,self._err])
			return False
		return True

	def add_argument(self, argument, value=None):
		argument = self._join(argument, value)
		self.command = self._join(self.command, argument)
		return

	def submit(self):
		self._run()
		return self._check_error()
	
	def output(self):
		self._run()
		return self._out if self._check_error() else None


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
### HELPER FUNCTIONS
################################################################################
def find( file, start='./' ):
	start += '/'
	if exists(start+file):
		return start+file
	for dir in filter(isdir, glob(start+'*')):
		return find(file, start=dir)
	return None


def get_rosetta_exe( exe, tools=False ):
	rosetta = expandvars('$ROSETTA')
	if '$' in rosetta:
		rosetta = '~/src/rosetta/'
	if tools:
		rosetta += '/tools/'
	else:
		rosetta += '/main/source/bin/'
	return find( exe, start=rosetta )


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


def get_score_data( filename, colnames=['score'], sort=None, filters=None, keep=None ):
	if not isinstance(colnames, list):
		colnames = [ colnames ]
	data = []
	colidx = None
	with open( filename, 'r' ) as f:
		for line in f:
			if not "SCORE:" in line:
				continue
			cols = line.split()
			if "description" in line:
				colidx = map(cols.index, filter(cols.count, colnames))
				continue
			data.append(tuple([float(cols[idx]) for idx in colidx]))
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


def get_working_target():
	return basename(os.getcwd())


def get_silent_file( filename=['region_FINAL.out','swm_rebuild.out'], dir=None ):
	if not isinstance(filename, list):
		filename = [ filename ]
	if dir is not None:
		filename = ['/'.join([dir,file]) for file in filename]
	silent_files = filter(exists, filename)
	return silent_files[0] if len(silent_files) else None


def get_native_pdb():
	target = get_working_target()
	native_pdb = glob( target+'_????_RNA.pdb' )[0]
	return native_pdb


def get_start_pdb_list():
	target = get_working_target()
	start_pdbs = glob( target+'_START*_????_RNA.pdb' )
	return start_pdbs


def get_motif_length():
	length = len(''.join(get_sequences(get_native_pdb())[0]))
	for start_pdb in get_start_pdb_list():
		length -= len(''.join(get_sequences(start_pdb)[0]))
	return length


def get_pdb_id():
	return get_native_pdb().split('_')[-2].upper()


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


def get_flag( flag ):
	if not exists( 'flags' ):
		return None
	with open( 'flags', 'r' ) as f:
		for line in f:
			if flag not in line: 
				continue
			return line.strip()
	return None 	


def virtualize_missing_residues( silent_file ):
	silent_file_out = silent_file.replace(".out","_full_model.out")
	if exists( silent_file_out ):
		Command( "rm -f ", args=silent_file_out ).submit()
	build_full_model_exe = get_rosetta_exe( "build_full_model" )
	weights = get_flag( "-score:weights" ).split(' ')[-1]
	command = Command( build_full_model_exe )
	command.add_argument( "-in:file:silent", value=silent_file )
	command.add_argument( "-out:file:silent", value=silent_file_out )
	if weights is not None:
		command.add_argument( "-score:weights", value=weights )
	command.add_argument( "-virtualize_built", value="true" )
	command.submit()
	return silent_file_out


def create_cluster_silent_file( silent_file ):
	if 'swm' in silent_file:
		silent_file = virtualize_missing_residues( silent_file )
	cluster_rmsd = 2.0 
	suite_cluster_rmsd = 2.5 
	no_graphic = False
	ignore_unmatched_virtual_res = False
	common_args_file = None 
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
	if not no_graphic:
		command.add_argument( "-no_graphic", value="False" )
	if ignore_unmatched_virtual_res:
		command.add_argument( "-ignore_unmatched_virtual_res", value="True" )
	if common_args_file:
		command.add_argument( "-common_args", value=common_args_file )
	command.add_argument( " > LOG_create_cluster_silent_file.txt" )
	command.submit()	
	return cluster_silent_file


def get_lowest_energy_cluster_centers( nclusters=5 ):
	silent_file = get_silent_file()	
	cluster_silent_file = create_cluster_silent_file( silent_file )
	if not exists( cluster_silent_file ):
		print "CLUSTER_SILENT_FILE NOT FOUND FOR TARGET: %s" % target
		print "USING SILENT_FILE: %s/%s" % (target, silent_file)
		cluster_silent_file = silent_file
	cluster_center_list = []
	score_types = ['score', get_rmsd_type(cluster_silent_file)]
	data = get_score_data( cluster_silent_file, colnames=score_types, sort='score', keep=nclusters )
	if data is None:
		return cluster_center_list
	for idx, (energy, rmsd) in enumerate(data, start=1):
		if idx > nclusters:
			break
		cluster_center_list.append( [idx, rmsd, energy] )
	return cluster_center_list


################################################################################
### MAIN FUNCTIONS
################################################################################
def get_target_properties():
	''' 
		column:	| Motif Properties |
		subcolumn: | Length | PDB	 |

	'''
	length = get_motif_length()
	pdb_id = get_pdb_id()
	return [ length, pdb_id ]


def get_best_of_lowest_energy_cluster_centers():
	'''
		column:	| Best of Five Lowest Energy Cluster Centers |
		subcolumn: | Cluster Rank | RMSD | Rosetta Energy (RU)  |

	'''
	cluster_centers = get_lowest_energy_cluster_centers()
	for idx, cluster_center in enumerate(cluster_centers):
		rmsd = cluster_center[1]
		if idx and rmsd >= best_rmsd:
			continue
		cluster_centers[0] = cluster_center
		best_rmsd = rmsd
	return cluster_centers[0]


def get_lowest_rmsd_model():
	'''
		column:	| Lowest RMSD Model |
		subcolumn: | RMSD			  |

	'''
	silent_file = get_silent_file()
	rmsd_type = get_rmsd_type( silent_file )
	rmsd = get_score_data( silent_file, colnames=rmsd_type, sort=rmsd_type, keep=1 )
	return [ rmsd ]


def get_lowest_energy_sampled( opt_exp_inpaths ):
	'''
		column:	| Lowest Energy Sampled						 |
		subcolumn: | Rosetta Energy (RU) | E-Gap to Opt. Exp. (RU) |

	'''
	silent_file = get_silent_file()
	energy = get_score_data( silent_file, sort='score', keep=1 )
	if energy is None:
		return [ None, None ]
	opt_exp_energy = get_opt_exp_score( opt_exp_inpaths )
	if opt_exp_energy is None:
		return [ energy, None ]
	energy_gap = energy - opt_exp_energy
	return [ energy, energy_gap ]
