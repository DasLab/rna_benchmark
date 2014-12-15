#!/usr/env/python

from os import popen
command = "grep 'seconds' KEEP_LOG_FILE/CONDOR/*/*/*/*.out"

lines = popen( command ).readlines()

total_time = 0
for line in lines:
	cols = line.split()
	for col in cols:
		if '.' not in col: continue 
		try:
			time = float(col)
			total_time+=time
		except:
			continue
		#print col

print 'Total Time: ', total_time, ' seconds ' 
print 'Total Time: ', (total_time/60.), ' minutes ' 
print 'Total Time: ', ((total_time/60.)/60.), ' hours ' 
print 