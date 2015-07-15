#!/usr/bin/python

import rnamake.motif_factory as motif_factory
import sys

#m = motif_factory.factory.motif_from_file("1gid_RNA_119-126_196-202.out.1.pdb")
input_pdb = sys.argv[1]
m = motif_factory.factory.motif_from_file(input_pdb)

print m.score
