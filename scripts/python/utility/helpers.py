#!/usr/bin/python

###############################################################################
### import modules
###############################################################################
import os
import sys
import subprocess

###############################################################################
### helper function
###############################################################################
def init_environ( options ):
    if 'rosetta' in dir(options) and os.path.exists(options.rosetta):
        os.environ['ROSETTA'] = options.rosetta
    if not 'ROSETTA' in os.environ or not os.path.exists(os.environ['ROSETTA']):
        print 'WARNING:', os.environ['ROSETTA'], 'is not a working Rosetta repository!!!'
        print 'The path to Rosetta must be defined in $ROSETTA, or specified via --rosetta'
        print 'To set this variable, add the following to your .bashrc or .zshrc:'
        print 'export ROSETTA=/path/to/rosetta/\n'
    assert('ROSETTA' in os.environ)
    assert(os.path.exists(os.environ['ROSETTA']))
    os.environ['ROSETTA_BIN'] = os.environ['ROSETTA'] + '/main/source/bin/'
    os.environ['ROSETTA_DB'] = os.environ['ROSETTA'] + '/main/database/'
    os.environ['ROSETTA_DB_WEIGHTS'] = os.environ['ROSETTA_DB'] + '/scoring/weights/'
    os.environ['SWA_DAGMAN_TOOLS'] = os.environ['ROSETTA'] + '/tools/SWA_RNA_python/SWA_dagman_python/'
    return True 

###############################################################################
def safe_submit( command, allow_retry=False, max_retry=3 ):
    if isinstance(command, str):
        command = command.split()
    assert( isinstance(command, list) )
    if not allow_retry: 
        max_retry = 1
    for attempt in xrange(max_retry):
        stdout, stderr = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        ).communicate()
        if not stderr or not len(stderr):
            return stdout
        print "STDOUT:", stdout
        print "STDERR:", stderr
    return -1

###############################################################################
def parse_flags(string_in, replacements = {}):
    if isinstance(string_in, list):
        string_in = ' '.join(string_in)
    for old, new in replacements.iteritems():
        string_in = string_in.replace(old, new)
    substrings = [s for s in string_in.split(' ') if s not in ['','-','#']]
    flags = []
    for idx, substring in enumerate(substrings):
        if substring.startswith('-'):
            flags.append([substring])
            continue
        if not len(flags):
            continue
        flags[-1].append(substring)
    flags = dict([(f.pop(0), ' '.join(f)) for f in flags])
    return flags
         
