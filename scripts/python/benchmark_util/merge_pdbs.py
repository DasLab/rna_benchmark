#!/usr/bin/python

###############################################################################
### imports
###############################################################################
from sys import argv, exit
from os.path import exists

###############################################################################
### helper functions
###############################################################################
def merge_filenames(filenames):
    ext = [f[f.index('.'):] for f in filenames][0]
    return '_'.join([f.replace(ext,'') for f in filenames])+ ext

###############################################################################
### main function
###############################################################################
def merge_pdbs(input_pdbs, overwrite = False):
    output_pdb = merge_filenames( input_pdbs )
    if exists(output_pdb) and overwrite is False:
        return output_pdb
    fid_out = open(output_pdb,'w')
    atomno, chain, resno = None, None, None
    prev_atomno, prev_chain, prev_resno = None, None, None
    for pdb_idx, input_pdb in enumerate(input_pdbs):
        input_lines = open(input_pdb,'r').readlines()
        atom_idx, res_idx = 0, 0
        for line_idx, line in enumerate(input_lines):
            line = line.strip()
            if 'ATOM' not in line:
                continue
            atomno, chain, resno = int(line[7:11]), line[21], int(line[22:26])
            if pdb_idx:
                atom_idx += 1
                if line_idx:
                    res_idx += ( resno != int(input_lines[line_idx-1][22:26]) )
                else:
                    res_idx += 1
                atomno = prev_atomno + atom_idx
                chain = prev_chain
                resno = prev_resno + res_idx
                fixed_line = line[:7]
                fixed_line += '%4d' % atomno
                fixed_line += line[11:21]
                fixed_line += '%s' % chain
                fixed_line += '%4d' % resno
                fixed_line += line[26:]
                line = fixed_line
            fid_out.write(line+'\n')
        prev_atomno = atomno
        prev_chain = chain
        prev_resno = resno
    fid_out.close()
    return output_pdb

###############################################################################
### main
###############################################################################
if __name__=='__main__':
    
    overwrite = any ( '-overwrite' in arg for arg in argv )
    if overwrite:
        argv.remove('-overwrite')

    if len(argv) < 3:
        print "Must provide atleast 2 existing PDBs"
        exit(1)

    input_pdbs = argv[1:]
    assert( all ( exists(f) for f in input_pdbs ) )
    
    print merge_pdbs( input_pdbs, overwrite )



