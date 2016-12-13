#!/usr/bin/python
import rnamake
from parse_tag import parse_tag
from make_tag import make_tag_with_dashes_and_commas
import rnamake.motif_type as motif_type
import sys, os
from get_surrounding_res import get_surrounding_res_tag
from get_sequence import get_sequences_for_res
from utility import file_handlers, info_handlers
import rnamake.pose_factory as pf

RNAMAKE=os.path.expandvars("$RNAMAKE")
assert(os.path.exists(RNAMAKE))



def usage():
    print "usage:", __file__, " [-radius <RADIUS>] [-design_bps]"

try:
    input_pdb = sys.argv[1]
except:
    usage()

design_bps = "-design_bps" in sys.argv

radius = 15
if '-radius' in sys.argv:
    radius_idx = sys.argv.index('-radius')
    radius = int(sys.argv[radius_idx+1])
    print "using radius = ",radius



info_fid = file_handlers.TargetDefinitionsFile()

#Read in pose
#p = rnamake.pose.Pose(RNAMAKE+"/examples/getting_started/resources/p4p6")
#p = rnamake.pose.Pose(pdb=input_pdb)
p = pf.factory.pose_from_file(input_pdb)

#loop over motifs
for iterator,motif in enumerate(p.motifs(motif_type.ALL),start=1):

    if motif.mtype == motif_type.HELIX:
        continue
    
    
    # figure out target attributes
    motif_res = []
    
    name = input_pdb.replace(".pdb", "")

    for chains in motif.chains():
        first = chains.first().num
        last = chains.last().num
        
        if not design_bps:
            first += 1
            last -= 1
            
        restag = "{}".format(first)
        if first != last:
            restag += "-{}".format(last)
    
        name += "_{}".format(restag)
        motif_res.append('{}:{}'.format(chains.first().chain_id,restag))  

    print name

    motif_res = ','.join(motif_res)
    print "motif_res = ", motif_res
    
    input_res = get_surrounding_res_tag(input_pdb, motif_res, radius, csv=True, verbose=True)
    
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
    td.name = name
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
