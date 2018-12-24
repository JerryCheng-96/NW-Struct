#!/usr/bin/python3
import re
import sys
import pickle
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP


print("Dump_DSSP.py - Dump the parsed DSSP result list to binary file using Library pickle.")
if len(sys.argv) < 2:
    print("Usage: Dump_DSSP.py input_filename [output_filename] [struct_no]")
    exit(0)

pp = PDBParser(PERMISSIVE=1)
filename = sys.argv[1]
try:
    pdb_id = re.findall(r"([^/]+?)\.pdb$", filename)[0]
except IndexError as ie:
    print("The input file should be a .pdb.")
    raise(ie)
    
output_filename = pdb_id + "_DSSP.dat"
struct_no = 0 
if len(sys.argv) >= 3:
    output_filename = sys.argv[2]
if len(sys.argv) >= 4:
    struct_no = sys.argv[3]

struct = pp.get_structure(pdb_id, filename)

try:
    model = struct[struct_no]
except IndexError as ie:
    print("No such structure with number \"" + struct_no + "\"")
    raise(ie)

dssp = DSSP(model, filename, dssp="mkdssp")

f = open(output_filename, "wb")
pickle.dump(dssp, f)
print("Dumped DSSP result list for PDB file " + filename + " to " + output_filename)