#!/usr/bin/python

###############################################################################
### import modules
###############################################################################
import sys
import os
import info_handlers

###############################################################################
### CLASSES
###############################################################################
class TargetDefinitionsFile(object):
    
    def __init__(self):
        self.target_definitions = []

    def add_target_definition(self, target_definition):
        assert( target_definition.name != td.name for td in self.target_definitions )
        self.target_definitions.append(target_definition)
        return True

    def load(self, fid):
        attrs, values = None, None
        if not isinstance(fid, file):
            fid = open(fid)
        for line in fid:
            # parse/check values
            values = line.split('#')[0].strip().split('\t')
            if len(values) < 2:
                continue
            # header line
            if attrs is None:
                attrs = [v.lower() for v in values]
                continue
            # new target definition
            target_definition = info_handlers.TargetDefinition()
            for attr, value in zip(attrs, values):
                if attr in 'sequence':
                    value = value.lower()
                setattr(target_definition, attr, value)
            target_definition.finalize()
            print "Loading TargetDefinition for:", target_definition.name
            self.add_target_definition( target_definition )
        fid.close()
        return True

    def _header(self):
        attrs = ['Name', 'Sequence', 'Secstruct', 'Working_res', 'Native', 'Input_res', 'Extra_flags']
        return '\t'.join(attrs)
        
    def save(self, fid):
        if not self.validate():
            return False
        if not isinstance(fid, file):
            fid = open(fid, 'w')
        fid.write(self._header() + '\n')
        fid.write('\n'.join([td._to_str(sep='\t') for td in self.target_definitions]))
        fid.close()
        return True

    def validate(self, verbose = False):
        if verbose is True:
            print '\n\n'.join([td._to_str(sep='\n') for td in self.target_definitions])
        for td in self.target_definitions:
            seqblocks = td.sequence.split(',')
            ssblocks = td.secstruct.split(',')
            wrblocks = td.working_res.split(',')
            if td.name is False:
                return False
            if not (len(seqblocks) == len(ssblocks) == len(wrblocks)):
                return False
            if td.native is None:
                return False
        return True

