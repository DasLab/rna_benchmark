#!/usr/bin/python
import sys, re, glob

n_cols=3
i=1
with open("README.md", "w") as out:
	out.write("<table>\n")
	for i, native_pdb in enumerate(glob.glob("*NATIVE*.pdb")):
		#native_pdb = re.sub("\_", "\\\_", native_pdb)
		data_str = "<td text-align=\"center\">%s<br /><img src=%s /></td>" % (re.sub("_NATIVE.*", "", native_pdb), re.sub(".pdb", ".png", native_pdb))
		if i % 3 == 0:
			out.write( "<tr>%s" % data_str )
		elif i % 3 == 1:
			out.write( data_str )
		else:
			out.write( "%s</tr>\n" % data_str )
	out.write("</table>\n")	

        #echo "![native structure of ${NATIVE_PDB/NATIVE*//}](${NATIVE_PDB/.pdb/.png})" >> README.md
