#!/usr/bin/python
import rnamake
import rnamake.motif_type as motif_type
import sys, os
from utility import file_handlers, info_handlers

RNAMAKE=os.path.expandvars("$RNAMAKE")
assert(os.path.exists(RNAMAKE))

info_fid = file_handlers.TargetDefinitionsFile()
input_pdb = sys.argv[1]

#Read in pose
p = rnamake.pose.Pose(RNAMAKE+"/examples/getting_started/resources/p4p6")
#p=rnamake.pose.Pose(pdb=input_pdb)

#loop over motifs
for iterator,motif in enumerate(p.motifs.all_motifs,start=1):

    if motif.mtype == motif_type.HELIX:
        continue
    
    # figure out target attributes
    working_res = []
    for residue in motif.chains():
        working_res.append('%s:%d-%d'%(residue.first().chain_id,residue.first().num,residue.last().num))  
    working_res = ','.join(working_res)
    
    #construct target object
    td = info_handlers.TargetDefinition()
   
    # set target attributes
    td.name = input_pdb.replace(".pdb", '_'+str(iterator))
    td.sequence = motif.sequence().replace('&',',').replace('+',',')
    td.secstruct = motif.secondary_structure().replace('&',',').replace('+',',')
    td.native = input_pdb
    td.working_res = working_res

    info_fid.add_target_definition(td)

#save '2r8s.txt'
motif_file = os.path.basename(os.path.dirname(os.path.abspath(input_pdb)))+'.txt'
print 'Saving TargetDefinitions file:', motif_file

#write save method
#info_fid.save(open(motif_file,"w"))

print 'Validating', motif_file
assert( info_fid.validate(verbose = True) )
