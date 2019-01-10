import PDB_Parse as pdbPar
import C_Interface as c
import Scores as scr

from Bio.PDB.PDBIO import PDBIO

#struct_1 = pdbPar.Get_PDBStruct("pdbs/tm_model.pdb")[0]
#struct_2 = pdbPar.Get_PDBStruct("pdbs/tm_exp.pdb")[0]
#resis_1 = pdbPar.Get_AminoResidues(struct_1, chain_id=' ')
#resis_2 = pdbPar.Get_AminoResidues(struct_2, chain_id=' ')
#
#rotMat, transVtr, a, b = scr.TMscoreAligned(resis_1, resis_2, len(resis_1))
#for atom in struct_1.get_atoms():
#    atom.transform(rotMat.T, transVtr)
#io = PDBIO()
#io.set_structure(struct_1)
#io.save("rotated.pdb")


pAmino_1, seqStr_1, len_1 = c.PDB2CAminoPtrArray(pdbPar.Get_AnimoNVtr(pdbPar.Get_AminoResidues(pdbPar.Get_PDBStruct("pdbs/1MBA.pdb"))), pdbPar.Get_DSSP("pdbs/1MBA.pdb"))
pAmino_2, seqStr_2, len_2 = c.PDB2CAminoPtrArray(pdbPar.Get_AnimoNVtr(pdbPar.Get_AminoResidues(pdbPar.Get_PDBStruct("pdbs/101M.pdb"))), pdbPar.Get_DSSP("pdbs/101M.pdb"))
revTracePath = c.NW.Needleman_Wunsch(pAmino_1, pAmino_2, c.cos_simFunc, 10, 10, len_seq_1=len_1, len_seq_2=len_2)

struct_1 = pdbPar.Get_PDBStruct("pdbs/1MBA.pdb")
struct_2 = pdbPar.Get_PDBStruct("pdbs/101M.pdb")
resis_1 = pdbPar.Get_AminoResidues(struct_1, chain_id='A')
resis_2 = pdbPar.Get_AminoResidues(struct_2, chain_id='A')

resStr_1, resStr_2 = c.NW.TracePath2AlignedStr(revTracePath, seqStr_1.decode('ascii'), seqStr_2.decode('ascii'))
res_1, res_2 = c.NW.TracePath2OnlyAlignedList(revTracePath, resis_1, resis_2)

rotMat, transVtr, a, b = scr.TMscoreAligned(res_1, res_2, len_2)
for atom in struct_1[0].get_atoms():
    atom.transform(rotMat.T, transVtr)
io = PDBIO()
io.set_structure(struct_1[0])
io.save("rotated_1MBA.pdb")
print()
print(resStr_1)
print(resStr_2)
