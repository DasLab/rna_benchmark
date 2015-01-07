#!/usr/bin/python

from os import popen
from sys import argv

if len( argv ) > 1:
	targets = argv[1:]	
else:
	targets = ["."]

for target in targets:
	print '\n%s' % target
	for time_key in ['resources_used.cput', 'resources_used.walltime']:	
		command = "grep %s %s/LOG_DEBUG_QSTAT_FULL.txt" % (time_key,target)
		lines = popen( command ).readlines()
		total_seconds = 0.0
		for line in lines:
			if time_key not in line or '=' not in line: continue
			[hh,mm,ss] = line.split()[2].split(':')
			total_seconds += float(hh)*60.*60. + float(mm)*60. + float(ss)
		mm, ss = divmod(total_seconds, 60)
		hh, mm = divmod(mm, 60)
		time_string = "%d:%02d:%02d" % (hh, mm, ss)
		print '%-23s = %s' % (time_key, time_string)	


