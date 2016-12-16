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
        self.extra_flags_benchmark = {}
        self.input_res_benchmark = []
        self.extra_min_res_benchmark = []

    def add_target_definition(self, target_definition):
        assert( target_definition.name != td.name for td in self.target_definitions )
        self.target_definitions.append(target_definition)
        return True

    def load(self, fid):
        attrs, values = None, None
        if not isinstance(fid, file):
            fid = open(fid)
        # get all lines and iterate through -- this way we can check for
        # alternate formats...
        lines = fid.readlines()
        # 'old format' == characterized by 'Name\t' first characters
        if lines[0].startswith("Name\t"):
            for line in lines:
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
                    # AMW: Don't do this! Then we lose our carefully X[FOO]-ified sequence!
                    #if attr in 'sequence':
                    #    value = value.lower()
                    setattr(target_definition, attr, value)
                target_definition.finalize()
                print "Loading TargetDefinition for:", target_definition.name
                self.add_target_definition( target_definition )
        else:
            i = 0
            while i < len(lines):
                # Re-handle extra flags benchmark later
                if len(lines[i]) <= 4:
                    i += 1
                    continue
                elif lines[i].split()[0] == "Benchmark_flags:":
                    substrings = lines[i].split()[1:]
                    substrings = [s for s in substrings if s not in ['','-','#']]
                    flags = []
                    for idx, substring in enumerate(substrings):
                        if substring.startswith('-'):
                            flags.append([substring])
                            continue
                        if not len(flags):
                            continue
                        flags[-1].append(substring)
                    self.extra_flags_benchmark = dict([(f.pop(0), ' '.join(f)) for f in flags])
                    if '-input_res' in self.extra_flags_benchmark:
                        self.input_res_benchmark = extra_flags_benchmark.pop('-input_res')
                    if '-extra_min_res' in self.extra_flags_benchmark:
                        self.extra_min_res_benchmark = extra_flags_benchmark.pop('-extra_min_res')

                if lines[i].split()[0] == "Name:":
                    # new format
                    target_definition = info_handlers.TargetDefinition()
                    while i < len(lines) and len(lines[i]) > 4: # forgive a little extra whitespace
                        setattr(target_definition, lines[i].split()[0].split(':')[0].lower(), " ".join( lines[i].split()[1:] ) )
                        i += 1
                    target_definition.finalize()
                    print "Loading TargetDefinition for:", target_definition.name
                    self.add_target_definition( target_definition )

                i += 1
        fid.close()
        return True

    def _header(self):
        attrs = ['Name', 'Sequence', 'Secstruct', 'Working_res', 'Native', 'Input_res', 'Align_res', 'Extra_flags']
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
            if td.name is False:
                return False
            #seqblocks = td.sequence.split(',')
            #ssblocks = td.secstruct.split(',')
            #wrblocks = td.working_res.split(',')
            #if not (len(seqblocks) == len(ssblocks) == len(wrblocks)):
            #    return False
            if td.native is None:
                return False
        return True

