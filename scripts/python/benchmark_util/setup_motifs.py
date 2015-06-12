#!/usr/bin/python
import rnamake
from parse_tag import parse_tag
from make_tag import make_tag_with_dashes_and_commas
import rnamake.motif_type as motif_type
import sys, os
from get_surrounding_res import get_surrounding_res_tag
from get_sequence import get_sequences_for_res
from utility import file_handlers, info_handlers

RNAMAKE=os.path.expandvars("$RNAMAKE")
assert(os.path.exists(RNAMAKE))

info_fid = file_handlers.TargetDefinitionsFile()
input_pdb = sys.argv[1]

#Read in pose
p = rnamake.pose.Pose(RNAMAKE+"/examples/getting_started/resources/p4p6")
#p = rnamake.pose.Pose(pdb=input_pdb)

#loop over motifs
for iterator,motif in enumerate(p.motifs.all_motifs,start=1):

    if motif.mtype == motif_type.HELIX:
        continue
    
    print input_pdb.replace(".pdb", '_'+str(iterator))
    # figure out target attributes
    motif_res = []
    for residue in motif.chains():
        motif_res.append('%s:%d-%d'%(residue.first().chain_id,residue.first().num,residue.last().num))  
    
    motif_res = ','.join(motif_res)
    print "motif_res = ", motif_res
    
    input_res = get_surrounding_res_tag(input_pdb, motif_res, 15, csv=True, verbose=True)
    
    input_resnums, input_chains = parse_tag(input_res)
    motif_resnums, motif_chains = parse_tag(motif_res)
    print "input_res = ", input_res
    working_resnums = input_resnums+motif_resnums
    working_chains = input_chains+motif_chains

    working_res = [ (r,c) for r,c in zip(working_resnums, working_chains) ]
    working_res.sort(key=lambda x: x[0])

    working_resnums = [ r for r,c in working_res ]
    working_chains = [ c for r,c in working_res ]
    working_res = make_tag_with_dashes_and_commas(working_resnums,working_chains)
    
    print "working_res = ", working_res

    start_sequence = ''.join(get_sequences_for_res(input_pdb, working_res))
    sequences = []
    for block in working_res.split(','):
        block_res, block_chains = parse_tag(block)
        sequence = start_sequence[:len(block_res)]
        start_sequence = start_sequence[len(block_res):]
        sequences.append(sequence)
                                
    sequences=','.join(sequences)
    

    #construct target object
    td = info_handlers.TargetDefinition()
   
    # set target attributes
    td.name = input_pdb.replace(".pdb", '_'+str(iterator))
    #td.sequence = ','.join(get_sequences_for_res(input_pdb, working_res))
    td.sequence = sequences
    td.secstruct = ''.join([(x if x == ',' else '.') for x in td.sequence])
    td.native = input_pdb
    td.working_res = working_res
    td.input_res = input_res
    print td._to_str(sep='\n')
    info_fid.add_target_definition(td)

#save '2r8s.txt'
motif_file = '../'+os.path.basename(os.path.dirname(os.path.abspath(input_pdb)))+'.txt'
print 'Saving TargetDefinitions file:', motif_file

#write save method
info_fid.save(open(motif_file,"w"))

print 'Validating', motif_file
assert( info_fid.validate(verbose = True) )
