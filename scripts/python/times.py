#!/usr/bin/python

from os import popen
from sys import argv

targets = []
for arg in argv[1:]:
	targets.append( arg )

for target in targets:
	print '\nGetting times for %s' % target
	command = "grep 'seconds' %s/KEEP_LOG_FILE/CONDOR/*/*/*/*.out" % target
	lines = popen( command ).readlines()
	total_time = 0.0
	for line in lines:
		cols = line.split()
		for idx, col in enumerate(cols):
			if 'seconds' not in col:
				continue
			try:
				time = float(cols[idx-1])
			except:
				continue
			total_time+=time

	print 'Total Time: ', total_time, ' seconds '
	print 'Total Time: ', (total_time/60.), ' minutes '
	print 'Total Time: ', ((total_time/60.)/60.), ' hours '

