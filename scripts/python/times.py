#!/usr/bin/python

from os import popen
from sys import argv


if len( argv ) > 1:
	targets = argv[1:]
else:
	targets = ["."]

for target in targets:
	print '\n%s' % target
	
	# get number of loop and input residues
	for res_key in ['rmsd_res','input_res']:
		lines = open(target+'/README_SETUP.py', 'r').readlines()
		for line in lines:
			if res_key not in line: continue
			res_list = line.split(res_key)[1].replace("'","").replace('\n','').split()
			res_count = len(res_list)
			print '%-23s = %d' % (res_key, res_count)

	# get number of jobs submitted
	#njobs = len( open(target+'/ALL_SLAVE_JOB_IDS.txt', 'r' ).readlines() )
	#print '%-23s = %d' % ('num_slave_jobs', njobs)
	for sub_key in ['num_slave_nodes']:
		submit_args = open(target+'/SUBMIT_SWA', 'r').read().split('-')
		for arg in submit_args:
			if sub_key not in arg: continue
			sub_arg = arg.split()[1]
			print '%-23s = %s' % (sub_key, sub_arg)


	# get cpu and wall times
	for time_key in ['resources_used.cput', 'resources_used.walltime']:
		total_seconds = 0.0
		lines = open(target+'/LOG_DEBUG_QSTAT_FULL.txt', 'r').readlines()
		for line in lines:
			if time_key not in line or '=' not in line: continue
			[hh,mm,ss] = line.split()[2].split(':')
			total_seconds += float(hh)*60.*60. + float(mm)*60. + float(ss)
		if total_seconds == 0.0: continue
		mm, ss = divmod(total_seconds, 60)
		hh, mm = divmod(mm, 60)
		time_string = "%d:%02d:%02d" % (hh, mm, ss)
		print '%-23s = %s' % (time_key, time_string)


