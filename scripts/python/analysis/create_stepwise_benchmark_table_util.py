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
subcolumn_labels = {
	" "*10 : [
		"Motif Name"	
	],
	"Motif Properties" : [
		"Length",
		"PDB"
	],
	"Best of Five Lowest Energy Cluster Centers" : [
		"Cluster Rank",
		"RMSD",
		"Rosetta Energy (RU)"
	],
	"Lowest RMSD Model" : [
		"RMSD" 
	],
	"Lowest Energy Sampled" : [
		"Roestta Energy (RU)",
		"E-Gap to Opt. Exp. (RU)"
	]
}


################################################################################
### CLASSES
################################################################################
class Command(object):

	def __init__(self, command, args=None):
		self.command = self._lex(command, args)
		self.stdout = sp.PIPE
		self.stderr = sp.PIPE
		self.silent = False
		self._out = None
		self._err = None

	def _lex(self, command, args):
		if isinstance(command, list):
			command = ' '.join(command)
		if isinstance(args, list):
			args = ' '.join(args)
		return ' '.join([command, str(args)]) if args else command 

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
		argument = self._lex(argument, value) if value else argument
		self.command = self._lex(self.command, argument)
		return

	def submit(self):
		self._run()
		return self._check_error()
	
	def output(self):
		self._run()
		return self._out if not self._check_error() else None


class Table(object):

	def __init__(self, filename):
		self.filename = filename
		self.header = []
		self.column_labels = []
		self.subcolumn_labels = []
		self.rows = []
		self.delimiter = ' | '
		self.subdelimiter = '  '
		self.column_widths = []
		self.subcolumn_widths = []
		self.format_templates = []

	def save(self):
		with open( self.filename, 'w' ) as f:
			f.write( self.table_to_string() )
		return 

	def set_widths(self):
		self.set_column_widths()
		self.set_subcolumn_widths()
		return

	def set_column_widths(self):
		widths = [int(max(len(c),4)) for c in self.column_labels]
		self.column_widths = widths
		return

	def get_subcolumn_group_width(self, index):
		group_width = 0
		for n, label in enumerate(self.subcolumn_labels):
			if index != self.subcolumn_index[n]:
				continue
			group_width += len(label)
		return group_width

	def set_subcolumn_widths(self):
		if not len(self.subcolumn_labels): return
		widths = []
		for idx, label in zip(self.subcolumn_index, self.subcolumn_labels):
			weight = len(label) / float(self.get_subcolumn_group_width(idx))
			sub_width = np.floor( self.column_widths[idx] * weight )
			widths.append( sub_width )
		widths[0] = max([len(r[0]) for r in self.rows])
		widths = [int(max(w,4)) for w in widths]
		self.column_widths[0] = widths[0]
		self.subcolumn_widths = widths
		return 

	def to_string(self, value, width):
		if isinstance(value, int):
			return "%-*d" % (width, value)
		if isinstance(value, float):
			return "%-*.1f" % (width, value)
		value = str(value) if value else "--"
		return "%-*s" % (width, str(value))
		
	def row_to_string(self, row, widths, newline=True):
		column_strings = [self.to_string(c,w) for c,w in zip(row,widths)]
		if len(row) != len(self.column_labels):
			column_group = [[] for c in self.column_labels]
			for n, idx in enumerate(self.subcolumn_index):
				column_group[idx].append(column_strings[n])
			column_strings = [self.subdelimiter.join(g) for g in column_group]
		row_string = self.delimiter.join([s for s in column_strings]) 
		if newline:
			row_string += '\n'
		return row_string
	
	def table_to_string(self):
		self.set_widths()
		widths, subwidths = self.column_widths, self.subcolumn_widths
		table_string =  self.row_to_string(self.column_labels, widths)
		table_string += self.row_to_string(self.subcolumn_labels, subwidths)
		for row in self.rows:
			table_string += self.row_to_string(row, subwidths)
		return table_string

	def add_row(self, row):
		if not isinstance(row, list):
			row = [ row ]
		self.rows.append( row )
		return

	def add_column_labels(self, labels):
		self.column_labels = [l for l in labels]
		return

	def add_subcolumn_labels(self, labels):
		self.subcolumn_labels = []
		self.subcolumn_index = []
		for idx, colname in enumerate(self.column_labels):
			self.subcolumn_labels += labels[colname]
			self.subcolumn_index += [idx for x in labels[colname]]
		return

	def column_averages(self):
		col_accums = []
		col_counts = []
		for row_idx, row in enumerate(self.rows):
			idx = 0
			data = row[1:]
			for col in data:
				if isinstance(col, str):
					increment = 0.0
				else:
					col = float( col )
					increment = 1.0
				if row_idx == 0:
					col_accums.append( col )
					col_counts.append( increment )
				else:
					col_accums[idx] += col
					col_counts[idx] += increment
				idx += 1
		col_aves = []
		for accum, count in zip(col_accums,col_counts):
			if not count:
				col_aves.append(None)
				continue
			col_aves.append(accum/count)
		return col_aves


class TableRow(object):

	def __init__(self):
		self.columns = []

	def add_columns(self, cols):
		if not isinstance(cols, list):
			cols = [ cols ]
		for col in cols:
			self.columns.append( col )
		return

	def get_columns(self):
		return self.columns
	

################################################################################
### HELPER FUNCTIONS
################################################################################
def find( file, dir='./' ):
	dir += '/'
	if exists(dir+file):
		return dir+file
	for subdir in filter(isdir, glob(dir+'*')):
		return find(file, dir=subdir)
	return None


def get_rosetta_exe( exe, tools=False ):
	rosetta = expandvars('$ROSETTA')
	if '$' in rosetta:
		rosetta = '~/src/rosetta/'
	if tools:
		rosetta += '/tools/'
	else:
		rosetta += '/main/source/bin/'
	return find( exe, dir=rosetta )


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


def get_score_data( filename, colnames=['score'], sort=None, keep=None ):
	if not isinstance(colnames, list):
		colnames = [ colnames ]
	data = []
	colidx = None
	with open( filename, 'r' ) as f:
		for line in f:
			if not "SCORE:" in line:
				continue
			cols = line.split()
			if not colidx:
				colidx = map(cols.index, filter(cols.count, colnames))
				continue
			data.append(tuple([float(cols[idx]) for idx in colidx]))
	if sort and len(data):
		sort = colnames.index(sort) if isinstance(sort, str) else sort-1
		data = sorted(data, key=operator.itemgetter(sort))
	if len(colidx) == 1:
		data = [d[0] for d in data] 
	if keep and len(data):
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


def get_silent_file( filename=['region_FINAL.out','swm_rebuild.out'] ):
	if not isinstance(filename, list):
		filename = glob(filename)
	silent_files = filter(exists, filename)
	assert( len(silent_files) )
	return silent_files[0]


def get_native_pdb():
	target = basename(os.getcwd())
	native_pdb = glob( target+'_????_RNA.pdb' )[0]
	return native_pdb


def get_start_pdb_list():
	target = basename(os.getcwd())
	start_pdbs = glob( target+'_START*_????_RNA.pdb' )
	return start_pdbs


def get_motif_length():
	length = len(''.join(get_sequences(get_native_pdb())[0]))
	for start_pdb in get_start_pdb_list():
		length -= len(''.join(get_sequences(start_pdb)[0]))
	return length


def get_pdb_id():
	return get_native_pdb().split('_')[-2].upper()


def create_cluster_silent_file( silent_file ):
	cluster_rmsd = 1.5 
	suite_cluster_rmsd = 2.5 
	no_graphic = False
	ignore_unmatched_virtual_res = True
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
	score_types = ['score', get_rmsd_type(cluster_silent_file)]
	data = get_score_data( cluster_silent_file, colnames=score_types, sort='score', keep=nclusters )
	cluster_center_list = []
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
		column:    | Motif Properties |
		subcolumn: | Length | PDB     |

	'''
	length = get_motif_length()
	pdb = get_pdb_id()
	return [ length, pdb ]


def get_best_of_lowest_energy_cluster_centers():
	'''
		column:    | Best of Five Lowest Energy Cluster Centers |
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
		column:    | Lowest RMSD Model |
		subcolumn: | RMSD              |

	'''
	silent_file = get_silent_file()
	rmsd_type = get_rmsd_type( silent_file )
	rmsd = get_score_data( silent_file, colnames=rmsd_type, sort=rmsd_type, keep=1 )
	return [ rmsd ]


def get_lowest_energy_sampled():
	'''
		column:    | Lowest Energy Sampled                         |
		subcolumn: | Rosetta Energy (RU) | E-Gap to Opt. Exp. (RU) |

	'''
	silent_file = get_silent_file()
	energy = get_score_data( silent_file, sort=1, keep=1 )
	opt_exp_silent_file = get_silent_file( filename='*OPT_EXP.out' )
	if not opt_exp_silent_file:
		return [ energy, None ]
	opt_exp_energy = get_score_data( opt_exp_silent_file, sort=1, keep=1 )
	energy_gap = energy - opt_exp_energy
	return [ energy, energy_gap ]

