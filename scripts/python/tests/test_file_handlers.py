#!/usr/bin/python

###############################################################################
from utility import file_handlers, info_handlers


###############################################################################
# create TargetDefinitionFile object
print 'Creating TargetDefinitionsFile object'
info_fid = file_handlers.TargetDefinitionsFile()

# load file
info_file = 'test_data/test_TargetDefinitionsFile.txt'
print 'Loading TargetDefinitions file:', info_file
info_fid.load(open(info_file))

# validate
print 'Validating', info_file
assert( info_fid.validate(verbose = True) )

