#!/usr/bin/python

###############################################################################
### import modules
###############################################################################
import sys
import os
import helpers

###############################################################################
### CLASSES
###############################################################################
class TargetDefinition(object):

    def __init__(self, name = None):
        self.name = name
        self.sequence = None
        self.secstruct = None
        self.working_res = None
        self.native = None
        self.input_res = None
        self.extra_flags = None

    def ordered_attrs(self):
        return sorted(dir(self))

    def finalize(self):
        self.extra_flags = helpers.parse_flags(self.extra_flags)
        return True

    def _to_str(self, sep='\n'):
        attrs = [self.name, self.sequence, self.working_res, self.native, self.input_res, self.extra_flags]
        return sep.join([('-' if attr is None else attr) for attr in attrs]) 
        

###############################################################################
class Target(TargetDefinition):

    def __init__(self, definition):
        
        # inherit definition data
        self.name = definition.name
        self.sequence = definition.sequence
        self.secstruct = definition.secstruct
        self.working_res = definition.working_res
        self.native = definition.native
        self.input_res = definition.input_res
        self.extra_flags = definition.extra_flags
        
        # new member data
        self.fasta = None
        self.resnums = None
        self.chains = None
        self.helix_files = None
        self.working_native = None
        self.input_pdbs = None
        self.terminal_res = None
        self.extra_min_res = None
        self.loop_res = None
        self.VDW_rep_screen_pdb = None
        self.VDW_rep_screen_info = None
        self.align_pdb = None


