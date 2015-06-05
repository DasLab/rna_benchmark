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
         
