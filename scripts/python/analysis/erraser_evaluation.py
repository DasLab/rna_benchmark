#!/usr/bin/python

################################################################################
### import modules
################################################################################
import argparse
import subprocess
import multiprocessing as mp
import os
from os.path import exists,basename,dirname,isdir,isfile,realpath,expandvars
from sys import exit
from glob import glob

################################################################################
### GLOBAL DEFINITIONS
################################################################################
INPUT_PDB_TAG = '_phenix.pdb'
OUTPUT_PDB_TAG = '_phenix_erraser.pdb'
REFERENCE_PDB_TAG = '_phenix_erraser_ref.pdb'

################################################################################
### HELPER CLASSES
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

	def is_empty(self):
		return bool(not len(self._columns))

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

	def save(self):
		with open( self.filename, 'w' ) as f:
			f.write( self._table_to_string() )
		return 


################################################################################
### HELPER FUNCTIONS
################################################################################
def check_rosetta( path_to_rosetta ):
	if ( not len( path_to_rosetta ) ):
		path_to_rosetta = expandvars( '$ROSETTA' )
		if ( not len( path_to_rosetta ) or not exists( path_to_rosetta ) ):
			print
			print '############################################################'
			print 'WARNING: $ROSETTA must be defined as the path to rosetta!!!'
			print 'Export this variable, or put the following in your .bashrc:'
			print 'export ROSETTA=/path/to/rosetta/'
			print '############################################################'
			print
			exit(0)
	assert( exists( path_to_rosetta ) )
	return path_to_rosetta


def check_erraser_tools( path_to_rosetta ):
	path_to_rosetta = check_rosetta( path_to_rosetta )
	path_to_erraser = path_to_rosetta+'/tools/ERRASER/'
	assert( exists( path_to_erraser ) )

	origin_dir = os.getcwd()
	os.chdir( path_to_erraser )
	
	command = 'phenix.rna_validate'
	out, err = subprocess.Popen( command, shell=True,
								 stdout=subprocess.PIPE,
								 stderr=subprocess.PIPE ).communicate()
	if err and len(err):
		if not ('Usage:' in err or 'Usage:' in out):
			print 'COMMAND: '+command
			print 'OUTPUT: '+out
			print 'ERROR: '+err
			print
			print '############################################################'
			print 'Make sure you have downloaded and installed PHENIX!!!'
			print '1. Download PHENIX from http://www.phenix-online.org/.'
			print '2. Install PHENIX.' 
			print 
			print 'To check that you have properly installed PHENIX, run:' 
			print '$ phenix.rna_validate'
			print '############################################################'
			print 
			exit(0)
	
	command = './convert_to_phenix_python'
	out, err = subprocess.Popen( command, shell=True,
								 stdout=subprocess.PIPE,
								 stderr=subprocess.PIPE ).communicate()
	if (err and len(err)): 
		print 'COMMAND: '+command
		print 'OUTPUT: '+out
		print 'ERROR: '+err
		exit(0)

	os.chdir(origin_dir)
	return path_to_erraser


def warn(*args):
	print "[WARNING]", ' '.join(map(str, args))


################################################################################
### MAIN FUNCTION
################################################################################
def analyze(work_dir, (path_to_erraser)):

	############################################################################
	### change into directory
	############################################################################
	origin_dir = os.getcwd()
	if not isdir( work_dir ):
		warn( work_dir, "is not a directory" )
		return False
	os.chdir( work_dir )

	############################################################################
	### get exe and pdbs 
	############################################################################
	erraser_analysis_exe = path_to_erraser + '/erraser_analysis.py'
	try:
		input_pdb = glob( '????' + INPUT_PDB_TAG )[0]
		output_pdb = glob( '????' + OUTPUT_PDB_TAG )[0]
		reference_pdb = glob( '????' + REFERENCE_PDB_TAG )[0]
	except IndexError as e:
		warn( "PDBs not found for", work_dir )
		return False 

	############################################################################
	### run erraser_analysis on input and output pdbs
	############################################################################
	command = [erraser_analysis_exe, input_pdb, output_pdb]
	out, err = subprocess.Popen( command,
								 stdout=subprocess.PIPE, 
								 stderr=subprocess.PIPE
	).communicate()
	with open( 'erraser_analysis.out', 'w' ) as fid:
		fid.write( '%s %s %s' % (basename(command[0]),command[1],command[2]))
		fid.write( '\n\n' )
		fid.write( out )
	with open( 'erraser_analysis.err', 'w' ) as fid:
		fid.write( err )
	
	############################################################################
	### run erraser_analysis on input and reference pdbs
	############################################################################
	command = [erraser_analysis_exe, input_pdb, reference_pdb]
	out, err = subprocess.Popen( command, 
								 stdout=subprocess.PIPE, 
								 stderr=subprocess.PIPE 
	).communicate()
	with open( 'erraser_analysis_ref.out', 'w' ) as fid:
		fid.write( '%s %s %s' % (basename(command[0]),command[1],command[2]))
		fid.write( '\n\n' )
		fid.write( out )
	with open( 'erraser_analysis_ref.err', 'w' ) as fid:
		fid.write( err )
  

	############################################################################
	### write table info
	############################################################################
	results = {
		'input':[],
		'reference':[],
		'output':[]
	}
	for line in open('erraser_analysis_ref.out').readlines():
		cols = line.strip().split()
		if 'Clashscore' in line:
			results['input'].append(float(cols[1]))
			results['reference'].append(float(cols[2]))
			continue
		if 'Suite' in cols and 'Outlier' in cols:
			results['input'].append(int(cols[2]))
			results['reference'].append(int(cols[3]))
			continue
		if 'Pucker' in cols and 'Outlier' in cols:
			results['input'].append(int(cols[2]))
			results['reference'].append(int(cols[3]))
			continue
		if 'Bond' in cols and 'Outlier' in cols:
			results['input'].append(int(cols[2]))
			results['reference'].append(int(cols[3]))
			continue
		if 'Angle' in cols and 'Outlier' in cols:
			results['input'].append(int(cols[2]))
			results['reference'].append(int(cols[3]))
			break

	for line in open('erraser_analysis.out').readlines():
		cols = line.strip().split()
		if 'Clashscore' in line:
			results['output'].append(float(cols[2]))
			continue
		if 'Suite' in cols and 'Outlier' in cols:
			results['output'].append(int(cols[3]))
			continue
		if 'Pucker' in cols and 'Outlier' in cols:
			results['output'].append(int(cols[3]))
			continue
		if 'Bond' in cols and 'Outlier' in cols:
			results['output'].append(int(cols[3]))
			continue
		if 'Angle' in cols and 'Outlier' in cols:
			results['output'].append(int(cols[3]))
			break

	############################################################################
	### write table info
	############################################################################
	table_row = TableRow()
	table_row.add_columns( basename(work_dir) )
	table_row.add_columns( results['input'] )
	table_row.add_columns( results['reference'] )
	table_row.add_columns( results['output'] )


	############################################################################	
	### change back to working directory
	############################################################################
	os.chdir(origin_dir)
	return table_row


################################################################################
### main script
################################################################################
if __name__=='__main__':

	############################################################################
	### parse arguments
	############################################################################
	parser = argparse.ArgumentParser(
		description=(
			'Analyze ERRASER runs: check for improvement in clash score,' + 
			'rotameric conformations, and bond angle/distance outliers'
		)
	)

	parser.add_argument(
		'inpath',
		help='Path to run directory.'
	)
	parser.add_argument(
		'-t','--targets',
		nargs='+',
		help='List of targets.',
		default=[]
	)
	parser.add_argument(
		'--path_to_rosetta', 
		help='Path to working copy of rosetta.',
		default='',
	)
	parser.add_argument(
		'-j','--nproc',
		help='Number of jobs to run in parallel',
		default=(mp.cpu_count()-1)
	)
	parser.add_argument(
		'-f','--force',
		help='Rewrite tables even if they alread exist.',
		action='store_true'
	)
	options = parser.parse_args()

	inpath = options.inpath
	user_targets = sorted( options.targets )
	path_to_rosetta = options.path_to_rosetta
	nproc = int( options.nproc )

	########################################################################
	### checks and initializations 
	########################################################################
	assert( exists( inpath ) )
	origin_dir = os.getcwd()
	path_to_erraser = check_erraser_tools( path_to_rosetta )

	########################################################################
	### change into inpath
	########################################################################
	os.chdir( inpath )
	print
	print 'Running ERRASER analysis...'
	print
	print 'Inpath:'
	print inpath
	print 

	########################################################################
	### find all targets in inpath
	########################################################################
	targets = sorted(filter(lambda x: isdir(x), glob('*')))
	if len(user_targets):
		targets = filter(lambda x: x in user_targets, targets)
	print 'Targets:'
	for t in targets:
		print t
		print

	########################################################################
	### run erraser_analysis on all targets found in inpath
	########################################################################
	args = ( path_to_erraser )
	if nproc > 1:
		
		###################################################################
		### Parallel version
		###################################################################
		pool = mp.Pool(processes=nproc)
		
		out = [pool.apply_async(analyze, args=(t,args)) for t in targets]
		table_rows = [o.get() for o in out]
		
		pool.close()
		pool.join()
	   
	else:
		
		###################################################################
		### Serial Version
		###################################################################
	   	table_rows = [analyze(t,args) for t in targets]		   

	########################################################################
	### change back into working directory
	########################################################################
	os.chdir( origin_dir )

	########################################################################
	### create table for inpath
	########################################################################
	table_name = inpath.upper() + '.tab'
	if exists( table_name ) and not options.force:
		print "Table:", table_name, "already exists!!!"
	else:
		column_labels = [
			'Input PDB',
			'ERRASER - Ref',
			'ERRASER - Out'
		]
		subcolumn_labels = ['Target']
		subcolumn_labels += [ 
			'Clashscore', 
			'Suite Outlier',
			'Pucker Outlier',
			'Bond Outlier', 
			'Angle Outlier'
		] * 3
		table = Table( table_name )
		table.add_row( inpath )
		table.add_row( column_labels )
		table.add_row( subcolumn_labels )
		for table_row in table_rows:
			table.add_data_row( table_row.columns() )
		table.save()	

	########################################################################
	### process or print output
	########################################################################
	print 
	print 'Status:'
	for idx, (t, tr) in enumerate(zip(targets,table_rows)):
		status = '-' if tr.is_empty() else '+'
		print '[%s] %s' % (status, t)
	