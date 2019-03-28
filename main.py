import C_Interface as c
import Scores as scr

import PDB_Helper as ph

pdbFile1 = ph.ReadPDBAsAtomsList("pdbs/4x9h.pdb")
pdbFile2 = ph.ReadPDBAsAtomsList("pdbs/4x9h_sec.pdb")
resis_1 = ph.GetChain(pdbFile1, 'A').residues
resis_2 = ph.GetChain(pdbFile2, 'A').residues


pAmino_1, seqStr_1, len_1 = c.PDB2CAminoPtrArray(ph.Get_AminoNVtr(ph.GetChain(pdbFile1, 'A')))
pAmino_2, seqStr_2, len_2 = c.PDB2CAminoPtrArray(ph.Get_AminoNVtr(ph.GetChain(pdbFile2, 'A')))
revTracePath = c.NW.Needleman_Wunsch(pAmino_1, pAmino_2, c.cos_simFunc, 50, 0, len_seq_1=len_1, len_seq_2=len_2)

resStr_1, resStr_2 = c.NW.TracePath2AlignedStr(revTracePath, seqStr_1.decode('ascii'), seqStr_2.decode('ascii'))
res_1, res_2 = c.NW.TracePath2OnlyAlignedList(revTracePath, resis_1, resis_2)

rotMat, transVtr, a, b = scr.TMscoreAligned(res_1, res_2, len_2)

theChain = ph.GetChain(pdbFile1, 'A')
transAtomList = []
for residue in theChain.residues:
    for atom in residue.atoms.values():
        transAtomList.append(atom.transform(rotMat.T, transVtr).infoList)

ph.SavePDBFile('rotated.pdb', [ph.GetChain(transAtomList, 'A')])

print()
print(resStr_1)
print(resStr_2)
