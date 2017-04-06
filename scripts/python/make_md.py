#!/usr/bin/python
import sys, re, glob

n_cols=3
i=1
with open("README.md", "w") as out:
	out.write("<table>\n")
	for i, native_pdb in enumerate(glob.glob("*NATIVE*.pdb")):
		native_pdb = re.sub("\_", "\\\_", native_pdb)

		if i % 3 == 0:
			out.write( "<tr><td>%s</td>" % re.sub("NATIVE.*", "", native_pdb))
		elif i % 3 == 1:
			out.write( "<td>%s</td>" % re.sub("NATIVE.*", "", native_pdb))
		else:
			out.write( "<td>%s</td></tr>\n" % re.sub("NATIVE.*", "", native_pdb))
	out.write("</table>\n")	

        #echo "![native structure of ${NATIVE_PDB/NATIVE*//}](${NATIVE_PDB/.pdb/.png})" >> README.md
