#!/usr/bin/python
from __future__ import print_function
import sys, re, glob, os

# first, get descriptions
pwd = os.getcwd()
print(pwd)

descriptions_for_puzzle = {}
with open("%s.txt" % pwd) as input_file:
	lines = input_file.readlines()
	# Let's just assume description comes after name. 
	# In fact, it should be penultimate.
	remember_key = None
	for i, line in enumerate(lines):
		if line[0:4] == "Name": remember_key = line.split(": ")[1].strip()
		if line[0:11] == "Description":
			try:
				descriptions_for_puzzle[remember_key] = line.split(": ")[1].strip()
				remember_key = None
			except:
				continue


n_cols=3
i=1
with open("README.md", "w") as out:
	out.write("<table>\n")
	for i, native_pdb in enumerate(glob.glob("*NATIVE*.pdb")):
		#native_pdb = re.sub("\_", "\\\_", native_pdb)
		puzzle_name = re.sub("_NATIVE.*", "", native_pdb)
		data_str = "\t\t<td align=\"center\">%s<br /><img src=\"%s\" /><br />%s</td>\n" % (puzzle_name, re.sub(".pdb", ".png", native_pdb), descriptions_for_puzzle[puzzle_name])
		if i % 3 == 0:
			out.write( "\t<tr>\n%s" % data_str )
		elif i % 3 == 1:
			out.write( data_str )
		else:
			out.write( "%s\t</tr>\n" % data_str )
	# finish off the row if the glob has a remainder
	if len(glob.glob("*NATIVE*.pdb")) % 3 == 1:
		out.write("\t\t<td></td>\n\t\t<td></td>\n\t</tr>\n")
	elif len(glob.glob("*NATIVE*.pdb")) % 3 == 2:
		out.write("\t\t<td></td>\n\t</tr>\n")
	out.write("</table>\n")	

        #echo "![native structure of ${NATIVE_PDB/NATIVE*//}](${NATIVE_PDB/.pdb/.png})" >> README.md
