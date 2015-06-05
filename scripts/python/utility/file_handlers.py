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
                setattr(target_definition, attr, value)
            target_definition.finalize()
            print "Loading TargetDefinition for:", target_definition.name
            self.add_target_definition( target_definition )
        fid.close()
        return True

    def validate(self, verbose = False):
        for td in self.target_definitions:
            for attr in td.ordered_attrs():
                if attr.startswith('_'):
                    continue
                if verbose is True:
                    print attr + ': ', getattr(td, attr)
            print
        return True

