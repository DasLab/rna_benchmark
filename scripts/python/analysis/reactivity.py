#!/usr/bin/python

from read_pdb import read_pdb
import argparse
import numpy as np


parser = argparse.ArgumentParser(description='Predict DMS reactivity from pdb.')
parser.add_argument('pdb', default=None)
args = parser.parse_args()


( coords, pdb_lines, sequence, chains, residues ) = read_pdb( args.pdb )

#print 'coords: ', coords
#print 'sequence: ', sequence
#print 'chains: ', chains
#print 'residues: ', residues

reactive_res_atom_types = { 
	'A' : [ 'N1' ]#,
	#'C' : [ 'N3' ]
}


reactive_coords = {}

for idx, res in enumerate(residues):
	chain = chains[idx]
	res_type = sequence[ chain ][ res ].replace(' ', '')

	if res_type not in reactive_res_atom_types.keys(): continue
	
	res_coords = coords[ chain ][ res ]	
	for atom_type, atom_coords in res_coords.iteritems():
		atom_type = atom_type.replace(' ','')
		if atom_type in reactive_res_atom_types[ res_type ]:
			reactive_coords[ res ] = {}
			reactive_coords[ res ][ atom_type ] = atom_coords


protected_atoms = {}
reactive_atoms = {}
protection_cutoff = 3.0

for res, atom in reactive_coords.iteritems():
	protected_atoms[ res ] = {}
	reactive_atoms[ res ] = {}

	#print res, atom
	for atom_type, atom_coords in atom.iteritems():
		protected_atoms[ res ][ atom_type ] = {}
		reactive_atoms[ res ][ atom_type ] = {}

		is_protected = False 
		#print atom_type, atom_coords
		### See if these atom coords are protected
		### if other atom is within 3 angstroms, add to protected list and break
		
		for idx, nbr_res in enumerate(residues):
			
			if nbr_res == res: continue 
			if nbr_res-1 == res: continue 
			if nbr_res+1 == res: continue

			protected_atoms[ res ][ atom_type ][ nbr_res ] = {}
			reactive_atoms[ res ][ atom_type ][ nbr_res ] = {}

			chain = chains[idx]
			nbr_res_coords = coords[ chain ][ nbr_res ]	
			for nbr_atom_type, nbr_atom_coords in nbr_res_coords.iteritems():
				if 'H' not in nbr_atom_type: continue
				atom_coords = np.array( atom_coords )
				nbr_atom_coords = np.array( nbr_atom_coords )
				dist2 = np.sum((atom_coords-nbr_atom_coords)**2)
				if dist2 <= protection_cutoff**2 and dist2 > 0.0:
					protected_atoms[ res ][ atom_type ][ nbr_res ][ nbr_atom_type ] = np.sqrt(dist2) 
				else:
					reactive_atoms[ res ][ atom_type ][ nbr_res ][ nbr_atom_type ] = np.sqrt(dist2) 

protected_res = {}
obs_res_list = [3,4,5,6,11,12,13,14] 
for res,atom in protected_atoms.iteritems():
	for atom_type, nbr_residues in atom.iteritems():
		for nbr_res, nbr_atoms in nbr_residues.iteritems():
			for nbr_atom_type, dist in nbr_atoms.iteritems():

				#if not (res == 5 and nbr_res == 12) and not (res == 13 and nbr_res == 4): continue
				
				if res not in obs_res_list and nbr_res not in obs_res_list: continue  
				print res, atom_type, nbr_res, nbr_atom_type, dist
				protected_res[res] = nbr_res

#print 'protected_atoms: ', protected_atoms
#print 'reactive_atoms: ', reactive_atoms

#print 'REACTIVE COORDS: ', reactive_coords 


###from matplotlib import pyplot as plt 


###fig = plt.figure()
###ax = fig.add_subplot( 1, 1, 1 )

##xvar_data = [ res for res in residues ]
###yvar_data = [ res for res in reversed(residues) ]
###ax.plot( xvar_data, yvar_data, marker='.', markersize=4, linestyle=' ' )
###plt.show()