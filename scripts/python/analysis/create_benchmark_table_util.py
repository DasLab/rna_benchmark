#!/usr/bin/python

################################################################################
### import modules
################################################################################
import argparse
import subprocess as sp
import multiprocessing as mp
import os
from os.path import exists, basename, dirname, isdir, isfile, realpath, abspath
from glob import glob
from titles import get_title
from get_sequence import get_sequences


################################################################################
### GLOBAL DEFINITIONS
################################################################################
benchmark_dir = abspath(__file__).split('benchmark/')[0] + 'benchmark/'
column_names = [
	"          ",
	"Motif Properties",
	"Best of Five Lowest Energy Cluster Centers",
	"Lowest RMSD Model",
	"Lowest Energy Sampled" 
]
subcolumn_names = {
	"          " : [
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
        self.rows = []

    def save(self):
        with open( self.filename, 'w' ) as f:
            f.write( self.table_to_string() )
        return 

    def get_column_sizes(self):
        sizes = [len(col) for col in self.rows[0]]
        for row in self.rows:
            for i, col in enumerate(row):
                sizes[i] = max(len(col),sizes[i])
        return sizes

    def row_to_string(self, row, delim='\t'):
        sizes = self.get_column_sizes()
        colstrings = ['%-*s'%(size,col) for (size,col) in zip(sizes,row)]
        return delim.join([cs for cs in colstrings])
    
    def table_to_string(self):
        table_string = ''
        for row in self.rows:
            table_string += self.row_to_string(row, delim=' | ')
            table_string += '\n'
        return table_string

    def add_row(self, row):
        if not isinstance(row, list):
            row = [ row ]
        self.rows.append( row )
        return


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


def get_target_properties( target ):
	native_pdb = glob( target+'_????_RNA.pdb' )[0]
	start_pdbs = glob( target+'_START*_????_RNA.pdb' )
	# get length
	length = len(''.join(get_sequences(native_pdb)[0]))
	for start_pdb in start_pdbs:
		length -= len(''.join(get_sequences(start_pdb)[0]))
	# get PDB
	pdb = native_pdb.split('_')[-2].upper()
	return '%6d   %4s' % (length, pdb)


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

	#-write_score_only True only work in recalculate RMSD mode...
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
	native_pdb = glob( target+'_????_RNA.pdb' )[0]
	silent_files = filter(exists, silent_file_name)
	assert( len( silent_files ) )
	silent_file = silent_files[0]	
	cluster_silent_file = create_cluster_silent_file(native_pdb, silent_file)
	if exists(cluster_silent_file):
		print "FOUND:", cluster_silent_file
	cluster_center_list = []
	nclusters = 5
	for idx in xrange(1, nclusters+1):
		rmsd = 1.0*idx
		energy = 1.0*idx
		cluster_center_info = "%12d   %4.1f   %4.1f" % (idx, rmsd, energy)
		cluster_center_list.append( cluster_center_info )
	return cluster_center_list


def get_best_of_lowest_energy_cluster_centers( target, silent_file_name ):
	lowest_rmsd = 16.0
	best_cluster_center = ''
	for cluster_center in get_lowest_energy_cluster_centers( target, silent_file_name ):
		rmsd = float(cluster_center.split()[1])
		if rmsd >= lowest_rmsd:
			continue
		best_cluster_center = cluster_center
		lowest_rmsd = rmsd
	return best_cluster_center


def get_sorted_score_list( filename, colname='score' ):
	scores = []
	with open( filename, 'r' ) as f:
		colidx = None
		for line in f:
			if not "SCORE:" in line:
				continue
			cols = line.split() 
			if colname in line:
				colidx = cols.index( colname )
				continue
			if not colidx:
				continue
			scores.append( float( cols[colidx] ) )
	scores.sort()
	return scores


def get_lowest_rmsd_model( silent_file_name ):
	silent_files = filter(exists, silent_file_name)
	assert( len( silent_files ) )
	silent_file = silent_files[0]
	# get lowest rmsd 
	rmsd_type = 'rms_fill'
	if 'region_FINAL.out' in silent_file:
		rmsd_type = 'NAT_rmsd' 
	rmsd_list = get_sorted_score_list( silent_file, colname=rmsd_type )
	rmsd = rmsd_list[0]
	return "%2.2f" % (rmsd)


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
		return "%4.1f   --" % (energy)
	opt_exp_silent_file = opt_exp_silent_files[0]
	opt_exp_energy_list = get_sorted_score_list( opt_exp_silent_file )
	opt_exp_energy = opt_exp_energy_list[0]
	energy_gap = energy - opt_exp_energy
	return "%4.1f   %4.1f" % (energy, energy_gap)

