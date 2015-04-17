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
from os.path import exists, basename, dirname, isdir, isfile, realpath, abspath
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
### FUNCTIONS
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


def get_score_data( filename, colname='score' ):
	data = {}
	idx = 0
	with open( filename, 'r' ) as f:
		colidx = None
		for line in f:
			if not "SCORE:" in line:
				continue
			cols = line.split()
			if not colidx:
				if colname in cols:
					colidx = cols.index(colname)
				continue
			data[idx] = float(cols[colidx])  
			idx += 1
	return data


def get_sorted_score_data( filename, colname='score' ):
	data = get_score_data( filename, colname=colname )
	return dict(sorted(data.items(), key=operator.itemgetter(1)))


def get_sorted_score_list( filename, colname='score' ):
	data = get_sorted_score_data( filename, colname=colname )
	data = sorted(data.values())
	return data


def get_rmsd_type( silent_file ):
	rmsd_type = 'rms_fill'
	if 'region_FINAL.out' in silent_file:
		rmsd_type = 'NAT_rmsd' 
	return rmsd_type


def get_target_properties( target ):
	native_pdb = glob( target+'_????_RNA.pdb' )[0]
	start_pdbs = glob( target+'_START*_????_RNA.pdb' )
	# get length
	length = len(''.join(get_sequences(native_pdb)[0]))
	for start_pdb in start_pdbs:
		length -= len(''.join(get_sequences(start_pdb)[0]))
	# get PDB
	pdb = native_pdb.split('_')[-2].upper()
	return [ length, pdb ]


def create_cluster_silent_file( native_pdb, silent_file ):
	SWA_TOOLS="/home/calebgeniesse/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/"
	TOP_ENERGY_CLUSTERS_FOLDER = "TOP_ENERGY_CLUSTERS/"
	cluster_silent_file = "%s/top_energy_clusters.out" % TOP_ENERGY_CLUSTERS_FOLDER
	COMMON_ARGS_FILE = glob( "COMMON_ARGS/*%s" % silent_file )[0]
	CLUSTER_RMSD = 1.5 
	SUITE_CLUSTER_RMSD = 2.5 
	NO_GRAPHIC = False
	IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG = False
	IGNORE_UNMATCHED_VIRTUAL_RES = True
	if exists( cluster_silent_file ):
		return cluster_silent_file
	command = "mkdir -p %s " % TOP_ENERGY_CLUSTERS_FOLDER
	out, err = sp.Popen( command, shell=True,
						stdout=sp.PIPE, stderr=sp.PIPE ).communicate()
	if err and len(err):
		print err
	cluster_top_energy_command="%s/misc/SWA_cluster.py " % (SWA_TOOLS)
	cluster_top_energy_command+="-num_pose_kept 100 -distinguish_pucker false -extract_pdb False " 
	cluster_top_energy_command+="-cluster_rmsd %s " % (CLUSTER_RMSD)
	cluster_top_energy_command+="-suite_cluster_rmsd %s " % (SUITE_CLUSTER_RMSD)
	cluster_top_energy_command+="-silent_file %s "% (silent_file)
	cluster_top_energy_command+="-native_pdb %s "%(native_pdb)
	cluster_top_energy_command+="-common_args %s "%(COMMON_ARGS_FILE)
	cluster_top_energy_command+="-output_filename %s " %(cluster_silent_file)	
	cluster_top_energy_command+="-no_graphic %s " %(NO_GRAPHIC)
	cluster_top_energy_command+="-full_length_loop_rmsd_clustering True "
	if(IGNORE_FARFAR_NO_AUTO_BULGE_PARENT_TAG):
		cluster_top_energy_command+="-ignore_FARFAR_no_auto_bulge_parent_tag True "
	if(IGNORE_UNMATCHED_VIRTUAL_RES):
		cluster_top_energy_command+="-ignore_unmatched_virtual_res True "
	cluster_top_energy_command+=" > LOG_create_cluster_silent_file.txt"	
	command = cluster_top_energy_command
	out, err = sp.Popen( command, shell=True,
						stdout=sp.PIPE, stderr=sp.PIPE ).communicate()
	if err and len(err):
		print err
	if not exists( cluster_silent_file ):
		print "%s does not exist!!!" % cluster_silent_file
		return ''
	return cluster_silent_file


def get_lowest_energy_cluster_centers( target, silent_file_name ):
	nclusters = 5
	cluster_center_list = []
	native_pdb = glob( target+'_????_RNA.pdb' )[0]
	silent_files = filter( exists, silent_file_name )
	assert( len( silent_files ) )
	silent_file = silent_files[0]	
	cluster_silent_file = silent_file
	#######################################################
	### TODO: Get SWA_cluster.py to run
	#create_cluster_silent_file( native_pdb, silent_file )
	#######################################################
	if not exists( cluster_silent_file ):
		print "CLUSTER_SILENT_FILE NOT FOUND FOR TARGET: %s" % target
		print "USING SILENT_FILE: %s/%s" % (target, silent_file)
		cluster_silent_file = silent_file
	energy_data = get_sorted_score_data( cluster_silent_file )
	rmsd_type = get_rmsd_type( cluster_silent_file )
	rmsd_data = get_score_data( cluster_silent_file, colname=rmsd_type )
	for idx, energy in energy_data.iteritems():
		if idx >= nclusters:
			break
		rmsd = rmsd_data[idx]
		cluster_center_info = [ idx+1, rmsd, energy ]
		cluster_center_list.append( cluster_center_info )
	return cluster_center_list


def get_best_of_lowest_energy_cluster_centers( target, silent_file_name ):
	lowest_rmsd = 16.0
	best_cluster_center = None
	cluster_centers = get_lowest_energy_cluster_centers( target, silent_file_name )
	for cluster_center_info in cluster_centers:
		if not best_cluster_center:
			best_cluster_center = cluster_center_info
			continue
		rmsd = cluster_center_info[1]
		rmsd = float(rmsd) if rmsd else None
		if not rmsd or rmsd >= lowest_rmsd:
			continue
		best_cluster_center = cluster_center_info
		lowest_rmsd = rmsd
	return best_cluster_center


def get_lowest_rmsd_model( silent_file_name ):
	silent_files = filter(exists, silent_file_name)
	assert( len( silent_files ) )
	silent_file = silent_files[0]
	# get lowest rmsd 
	rmsd_type = get_rmsd_type( silent_file )
	rmsd_list = get_sorted_score_list( silent_file, colname=rmsd_type )
	rmsd = rmsd_list[0]
	return [ rmsd ]


def get_lowest_energy_sampled( silent_file_name ):
	silent_files = filter(exists, silent_file_name)
	assert( len( silent_files ) )
	silent_file = silent_files[0]
	# get lowest energy 
	energy_list = get_sorted_score_list( silent_file )
	energy = energy_list[0]
	# get energy gap
	opt_exp_silent_files = glob( '*OPT_EXP.out' )
	if not len( opt_exp_silent_files ):
		return ( energy, None )
	opt_exp_silent_file = opt_exp_silent_files[0]
	opt_exp_energy_list = get_sorted_score_list( opt_exp_silent_file )
	opt_exp_energy = opt_exp_energy_list[0]
	energy_gap = energy - opt_exp_energy
	return [ energy, energy_gap ]

