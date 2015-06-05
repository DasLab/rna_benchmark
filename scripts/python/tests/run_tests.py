#!/usr/bin/python

import sys
import os
import glob
import subprocess

path = os.path.dirname(os.path.abspath(__file__))
if not len(path):
    path = '.'

tests = glob.glob(path + '/test_*.py')

for test in tests:
    if path not in test:
        test = os.path.join(path, test)
    print "Running test:", test
    subprocess.call(['python', test])
