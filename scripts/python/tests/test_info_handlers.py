#!/usr/bin/python

###############################################################################
from utility import info_handlers, file_handlers

###############################################################################
# construct TargetDefinition object
target_def = info_handlers.TargetDefinition( 'test' )

# create TargetDefinitionFile object
print 'Creating TargetDefinitionsFile object'
info_fid = file_handlers.TargetDefinitionsFile()

# add TargetDefinition object 
print 'Adding TargetDefinition object'
info_fid.add_target_definition( target_def )

# validate
print 'Validating'
assert( info_fid.validate(verbose = True) )
