#!/usr/bin/python
import sys, re, glob

n_cols=3
i=1
with open("README.md", "w") as out:
	for i, native_pdb in enumerate(glob.glob("*NATIVE*.pdb")):
		native_pdb = re.sub("\_", "\\\_", native_pdb)

		if i % 3 == 0:
			out.write( "%s|" % re.sub("NATIVE.*", "", native_pdb))
		elif i % 3 == 1:
			out.write( "%s" % re.sub("NATIVE.*", "", native_pdb))
		else:
			out.write( "|%s\n" % re.sub("NATIVE.*", "", native_pdb))
			out.write("---------------------------------------------\n")
				

        #echo "![native structure of ${NATIVE_PDB/NATIVE*//}](${NATIVE_PDB/.pdb/.png})" >> README.md
