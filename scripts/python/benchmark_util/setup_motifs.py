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
for i,m in enumerate(p.motifs.all_motifs,start=1):

    if m.mtype == motif_type.HELIX:
        continue
    
    # figure out target attributes
    print m.sequence(), m.secondary_structure()
    # replace '&' with ',' and '+' with ',' in sequence and secondary_structure 

    working_res = []
    for c in m.chains():
        #working_res.append('%s:%d-%d'%(c.first().chain_id,c.first().num,c.last.num))  
        print c.first().num, c.last().num,c.first().chain_id

    working_res = ','.join(working_res)
    
    #construct target object
    td = info_handlers.TargetDefinition()
   
    # set target attributes
    # ex.td.sequence, td.secstruct, td.working_res
    # figure out name = input_pdb.replace(.pdb, '_'+str(i))
    td.name = input_pdb.replace(".pdb", '_'+str(i))
    td.sequence = m.sequence()
    td.sectruct = m.secondary_structure()
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
