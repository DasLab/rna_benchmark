#!/usr/bin/python

###############################################################################
### imports
###############################################################################
from sys import argv, exit
import subprocess as sp
import os
import os.path
import shutil
import glob

###############################################################################
### helper functions
###############################################################################
def safe_submit(command):
    out, err = sp.Popen(
        command,
        shell=isinstance(command, str),
        stdout=sp.PIPE,
        stderr=sp.PIPE
    ).communicate()
    if err:
        print out
        print err
        return False
    return True
       

###############################################################################
### main functions
###############################################################################
def combine_qsub_files_MPI():
    workdir = os.getcwd()
    mpi_batch_file = 'MPI_ONEBATCH.job'
    batch_files = glob.glob('%s/*/*/%s' % (workdir, mpi_batch_file))
    if not len(batch_files):
        batch_files = glob.glob('%s/*/%s' % (workdir, mpi_batch_file))
        if not len(batch_files):
            return False
    master_batch_files = []
    for batch_id, batch_file in enumerate(batch_files, start=1):
        job_lines = open( batch_file, 'r').readlines()
        for job_idx, job in enumerate(job_lines, start=1):
            jobid = (batch_id * job_idx) - 1
            if jobid % 16:
                id = len(master_batch_files)
                master_batch_file = mpi_batch_file.replace('.job','_%d.job'%id)
                master_batch_files.append(master_batch_file)
                open(mpi_batch_file, 'w')
            with open(mpi_batch_file, 'a') as master_batch_fid:
                master_batch_fid.write( job.strip()+'\n' )
    with open('qsubMPI_ONEBATCH', 'w') as master_submit_fid:
        for idx, master_batch_file in enumerate(master_batch_files):       
            job_file = 'qsubMPI_ONEBATCH_%d.sbatch' % idx
            with open(job_file, 'w') as fid:
                fid.write('#!/bin/bash\n')
                fid.write('#SBATCH -J mpi_onebatch%d\n' % idx)
                fid.write('#SBATCH -o %j.out\n')
                fid.write('#SBATCH -p normal\n')
                fid.write('#SBATCH -t 4:00:00\n')
                fid.write('#SBATCH -n 16\n')
                fid.write('#SBATCH -N 1\n')
                fid.write('#SBATCH -A TG-MCB120152\n')
                fid.write('pp_jobsub.py %s' % master_batch_file) 
                fid.write(' -cluster_name stampede')
                fid.write(' -nodelist $SLURM_NODELIST')
                fid.write(' -job_cpus_per_node $SLURM_JOB_CPUS_PER_NODE\n')
            master_submit_fid.write('sbatch %s\n' % job_file)
    return True


###############################################################################
### main
###############################################################################
if __name__=='__main__':
    
    import argparse

    parser = argparse.ArgumentParser(
        description='HPC Cluster Utility Tools.'
    )
    parser.add_argument(
        '--combine_qsub_files_MPI',
        help='Combine MPI job submission files from subdirectories.',
        action='store_true'
    )
    args = parser.parse_args()
    
    if args.combine_qsub_files_MPI:
        print "combine_qsub_files_MPI:", combine_qsub_files_MPI()
    
