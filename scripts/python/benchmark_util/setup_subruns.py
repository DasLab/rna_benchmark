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
### global vars
###############################################################################
SUBMIT_FILES = ['qsubMINI','qsubMPI','sbatchMINI']

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

def expand_tags(tags):
    if not isinstance(tags, list):
        tags = [ tags ]
    expanded_tags = []
    for idx, tag in enumerate(tags):
        expanded_tags.append(tag)
        if '-' not in tag:
            continue
        if not any ( c.isdigit() for c in tag ):
            continue
        idx = [tag.index(c) for c in tag if c.isdigit()][0]
        prefix = tag[:idx]
        [start, stop] = tag.replace(prefix, '').split('-')
        if not (start.isdigit() and stop.isdigit()):
            continue
        expanded_tags.pop()
        expanded_tags += [prefix+str(m) for m in range(int(start),int(stop)+1)]
    return expanded_tags
        

###############################################################################
### main functions
###############################################################################
def combine_qsubMPI():
    mpi_batch_file = 'MPI_ONEBATCH.job'
    batch_files= filter(os.path.isdir,glob.glob('./*/*/%s'%mpi_batch_file))
    master_batch_files = []
    for batch_id, batch_file in enumerate(batch_files, start=1):
        job_lines = open( batch_file, 'r').readlines()
        for job_idx, job in enumerate(job_lines, start=1):
            jobid = (batch_id * job_id) - 1
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

def setup_subruns(tags):
    readme_setup_fid = open('README_SETUP_MASTER', 'w')
    submit_fids = dict([(f, open(f,'w')) for f in SUBMIT_FILES]) 
    for tag in tags:
        workdir = os.getcwd()
        extra_flags = open('extra_flags_benchmark.txt.tmpl','r').readlines()
        rundir = os.path.basename(workdir)+'_'+tag.replace(':','')
        print "Setting up subdirectory:",rundir
        if os.path.exists(rundir):
            safe_submit('rm -r '+rundir)
        os.mkdir(rundir)
        os.chdir(rundir)
        if os.getcwd() == workdir:
            continue
        with open('extra_flags_benchmark.txt','w') as fid:
            for flag in extra_flags:
                flag = flag.replace('${VAR}', tag)
                flag = flag.replace('${GLOBAL_VARS}',','.join(tags))
                fid.write( flag.strip()+'\n' )
        with open('README_SETUP','w') as fid:
            args = open('../README_SETUP.tmpl','r').read().split(' ')
            for arg in args:
                if arg.startswith('../'):
                    arg = '../' + arg
                fid.write(arg+' ')
        os.chdir(workdir)
        readme_setup_fid.write(
            'cd %s; source ./README_SETUP; cd %s\n' % (rundir, workdir)
        )
        for submit_file, submit_fid in submit_fids.iteritems():
            submit_fid.write(
                'cd %s; source ./%s; cd %s\n' % (rundir, submit_file, workdir)
            )
    readme_setup_fid.close()
    for file, fid in submit_fids.iteritems():
        fid.close()
    return True 


###############################################################################
### main
###############################################################################
if __name__=='__main__':
    
    if len(argv) < 2:
        exit(1)
    
    if '-combine_qsubMPI' in argv[1]:
        combine_qsubMPI()
        exit(1)
    
    subrun_tags = argv[1:]
    subrun_tags = expand_tags(subrun_tags)
    setup_subruns(subrun_tags)
