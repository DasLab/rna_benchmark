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


################################################################################
### MAIN FUNCTION
################################################################################
def analyze(work_dir, (path_to_erraser)):

	############################################################################
	### change into directory
	############################################################################
	origin_dir = os.getcwd()
	os.chdir( work_dir )	

	############################################################################
	### get exe and pdbs 
	############################################################################
	erraser_analysis_exe = path_to_erraser + '/erraser_analysis.py'
	input_pdb = glob( '????' + INPUT_PDB_TAG )[0]
	output_pdb = glob( '????' + OUTPUT_PDB_TAG )[0]
	reference_pdb = glob( '????' + REFERENCE_PDB_TAG )[0]

	############################################################################
	### run erraser_analysis on input and output pdbs
	############################################################################
	command = [erraser_analysis_exe, input_pdb, output_pdb]
	out, err = subprocess.Popen( command,
								 stdout=subprocess.PIPE, 
								 stderr=subprocess.PIPE ).communicate()
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
								 stderr=subprocess.PIPE ).communicate()
	with open( 'erraser_analysis_ref.out', 'w' ) as fid:
		fid.write( '%s %s %s' % (basename(command[0]),command[1],command[2]))
		fid.write( '\n\n' )
		fid.write( out )
	with open( 'erraser_analysis_ref.err', 'w' ) as fid:
		fid.write( err )
  
	############################################################################	
	### change back to working directory
	############################################################################
	os.chdir(origin_dir)
	return True


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

	args = parser.parse_args()

	inpath = args.inpath
	user_targets = sorted( args.targets )
	path_to_rosetta = args.path_to_rosetta
	nproc = int( args.nproc )

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

	########################################################################
	### find all targets in inpath
	########################################################################
	targets = sorted(filter(lambda x: isdir(x), glob('*')))
	if len(user_targets):
		targets = filter(lambda x: x in user_targets, targets)
	print
	print 'Targets:'
	for t in targets:
		print t

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
		out = [o.get() for o in out]
		
		pool.close()
		pool.join()
	   
	else:
		
		###################################################################
		### Serial Version
		###################################################################
	   	out = [analyze(t,args) for t in targets]		   

	########################################################################
	### process or print output
	########################################################################
	print 
	print 'Status:'
	for idx, t in enumerate(targets):
		status = '+' if out[idx] else '-'
		print '[%s] %s' % (status, t)
	
	########################################################################
	### change back into working directory
	########################################################################
	os.chdir( origin_dir )
