#!/usr/bin/python

###############################################################################
### import modules
###############################################################################
import string
from subprocess import Popen, PIPE
from os.path import exists, expandvars
from os import system, uname
from make_tag import make_tag_with_dashes
from parse_options import get_resnum_chain
import subprocess

###############################################################################
### helper function
###############################################################################
def safe_submit( command, allow_retry=False, max_retry=3 ):
    if isinstance(command, str):
        command = string.split(command)
    assert( isinstance(command, list) )
    if not allow_retry: 
        max_retry = 1
    for attempt in xrange(max_retry):
        stdout, stderr = Popen(command, stdout=PIPE, stderr=PIPE).communicate()
        if not stderr or not len(stderr):
            return stdout
        print "STDOUT:", stdout
        print "STDERR:", stderr
    return -1

###############################################################################
def slice_out( inpath_dir, prefix, pdb, res_string, excise=False ):
    starting_native = inpath_dir+'/'+pdb
    assert( exists( starting_native ) )
    slice_pdb = prefix + pdb
    if not exists( slice_pdb ):
        if excise:  command = 'pdbslice.py %s -excise %s %s ' % ( starting_native, res_string, prefix )
        else:       command = 'pdbslice.py %s -subset %s %s ' % ( starting_native, res_string, prefix )
        system( command )
    assert( exists( slice_pdb ) )
    return slice_pdb

###############################################################################
def get_fullmodel_number( reschain, resnums, chains):
    for m in range( len( resnums ) ):
        if ( resnums[m] == reschain[0] ) and ( reschain[1] == '' or chains[m] == reschain[1] ): return m+1
    return 0

###############################################################################
def get_align_res( screen_pdb, working_pdb, working_fixed_res ):
    from read_pdb import read_pdb
    screen_align_res = []
    working_align_res = []
    for pdb in [ screen_pdb, working_pdb ]: assert( exists( pdb ) )
    #( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdb )
    screen_pdb_info = read_pdb( screen_pdb )
    working_pdb_info = read_pdb( working_pdb )
    for n in xrange( len( screen_pdb_info[4] ) ):
        screen_chain = screen_pdb_info[3][n]
        screen_res = screen_pdb_info[4][n]
        #working_reschain = pdb2pose( working_pdb_info[3], working_pdb_info[4], screen_chain, screen_res )   
        working_reschain = ( 0, '' )
        for m in xrange( len( working_pdb_info[3] ) ):
            if ( working_pdb_info[4][m] == screen_res ) and ( working_pdb_info[3][m] == screen_chain ): 
                working_reschain = ( working_pdb_info[4][m], working_pdb_info[3][m] )
        if ( working_reschain[0] > 0 ):# and ( working_reschain[0] in working_fixed_res ):
            screen_align_res.append( n )
            working_align_res.append( get_fullmodel_number(working_reschain,working_pdb_info[4],working_pdb_info[3]) ) 
    if len( screen_align_res ): screen_align_res_tag = make_tag_with_dashes( screen_align_res )
    else:   screen_align_res_tag = '0-0'
    if len( working_align_res ):    working_align_res_tag = make_tag_with_dashes( working_align_res )
    else:   working_align_res_tag = '0-0'
    return ( screen_align_res_tag, working_align_res_tag )

###############################################################################
def make_rna_rosetta_ready( pdb_file, sequence, reassign_chainids=True, allowed_chains=string.ascii_uppercase ):
    make_rna_rosetta_ready_cmdline = [ 'make_rna_rosetta_ready.py', pdb_file ]
    if reassign_chainids:
        chainid_str = ''
        for chainidx, seq in enumerate(sequence.split(',')):
            chainid_str += ''.join([ allowed_chains[chainidx] for res in seq ])
        make_rna_rosetta_ready_cmdline += [ '-reassign_chainids', chainid_str ]
    out, err = subprocess.Popen( make_rna_rosetta_ready_cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    rna_rosetta_ready_native = out.split()[-1] #native[ name ].lower().replace('.pdb', '_RNA.pdb')
    mv_cmdline = [ 'mv', rna_rosetta_ready_native, pdb_file ]
    out, err = subprocess.Popen( mv_cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()

###############################################################################
def make_replacements(s, replacements = {}):
    if isinstance(s, list):
        return [make_replacements(ss, replacements) for ss in s]
    for old, new in replacements.iteritems():
        s = s.replace(old, new)
    return s

###############################################################################
def parse_flags(substrings, replacements = {}):
    if isinstance(substrings, list):
        substrings = ' '.join(substrings)
    substrings = make_replacements(substrings, replacements)
    substrings = substrings.split(' ')
    substrings = [s for s in substrings if s not in ['','-','#']]
    flags = []
    for idx, substring in enumerate(substrings):
        if substring.startswith('-'):
            flags.append([substring])
            continue
        if not len(flags):
            continue
        flags[-1].append(substring)
    flags = dict([(f.pop(0), ' '.join(f)) for f in flags])
    return flags

###############################################################################
def merge_pdbs(pdbs_in, pdb_out):                                       
    if not isinstance(pdbs_in, list):                                   
        pdbs_in = [ pdbs_in ]           
    pdbs_in = ' '.join(pdbs_in)
    system('cat %s > %s' % (pdbs_in, pdb_out))              
    return pdb_out             

###############################################################################
def init_submit_files():
    hostname = uname()[1]
    submit_files = ['qsubMINI']
    if 'stampede' in hostname:
        submit_files = ['qsubMPI']
    if 'comet' in hostname:
        submit_files = ['sbatchMINI','qsubMPI']
    if 'sherlock' in hostname or 'SHERLOCK' in expandvars("$COMPUTER_CLUSTER_NAME"):
        submit_files = ['sbatchMINI','qsubMPI']
    for submit_file in submit_files:
        fid = open(submit_file, 'w')
        fid.close()
    return submit_files

###############################################################################
### CLASSES
###############################################################################
class FullModelInfo(object):
    
    def __init__(self, name = None):
        self.name = name
        self.resnums = None
        self.chains = None
             
    def set_resnums(self, resnums):
        self.resnums = resnums
        
    def set_chains(self, chains):
        self.chains = chains
        
    def extract_resnum_chains(self, resnum_chain_tag):
        resnums, chains = [], []
        if not isinstance(resnum_chain_tag, list):
            resnum_chain_tag = resnum_chain_tag.replace(';',',').split(',')
        for tag in resnum_chain_tag:
            get_resnum_chain(tag, resnums, chains)
        return resnums, chains

    def conventional_tag_to_full(self, resnum_chain_tag):
        resnums, chains = self.extract_resnum_chains(resnum_chain_tag)
        return self.conventional_to_full(zip(resnums, chains))
        
    def conventional_to_full(self, resnum_chain):
        '''
        Input: (res, chain) or [(res, chain), (res, chain)]
        '''
        if isinstance(resnum_chain, list):
            return map(self.conventional_to_full, resnum_chain)
        return get_fullmodel_number(resnum_chain, self.resnums, self.chains)
