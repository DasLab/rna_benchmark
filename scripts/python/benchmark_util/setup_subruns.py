#!/usr/bin/python

###############################################################################
### imports
###############################################################################
from sys import argv, exit
import subprocess as sp
import os
import os.path
import shutil

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

    subrun_tags = argv[1:]
    subrun_tags = expand_tags(subrun_tags)
    setup_subruns(subrun_tags)
