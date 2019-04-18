import C_Interface_New as c
import Scores as scr

import PDB_Helper as ph

#pdbFile1 = ph.ReadPDBAsAtomsList("pdbs/4x9h.pdb")
#pdbFile2 = ph.ReadPDBAsAtomsList("pdbs/4x9h_sec.pdb")
#pdbFile1 = ph.ReadPDBAsAtomsList("pdbs/Legacy/101M.pdb")
#pdbFile2 = ph.ReadPDBAsAtomsList("pdbs/Legacy/1MBA.pdb")
pdbFile1 = ph.ReadPDBAsAtomsList("pdbs/1atz.pdb")
pdbFile2 = ph.ReadPDBAsAtomsList("pdbs/1auo.pdb")
chain_1 = ph.GetChain(pdbFile1, 'A')
chain_2 = ph.GetChain(pdbFile2, 'A')
resis_1 = chain_1.residues
resis_2 = chain_2.residues
aminoNVtr_1 = chain_1.Get_AminoNVtr()
aminoNVtr_2 = chain_2.Get_AminoNVtr()
neighborVtr1 = chain_1.Get_SurroundVectorSet(4, 12)
neighborVtr2 = chain_2.Get_SurroundVectorSet(4, 12)

def CombinedSimFunc(ptr_1, ptr_2):
    ptr_1 = c.cast(ptr_1, c.POINTER(c.CPyAmino))
    ptr_2 = c.cast(ptr_2, c.POINTER(c.CPyAmino))
    aaNum_2 = ptr_1.contents.num
    aaNum_1 = ptr_2.contents.num
    if aaNum_1 == 0:
        print('Now At: (' + str(aaNum_1) + ', ' + str(aaNum_2) + ')')

    try:
        localSimVal = ph.Calc_VtrCosine(aminoNVtr_1[aaNum_1], aminoNVtr_2[aaNum_2])
        return localSimVal
        globalSimVal = ph.Calc_SimVectorSet(neighborVtr1[aaNum_1], neighborVtr2[aaNum_2])
    except IndexError:
        print(str(aaNum_1) + ', ' + str(aaNum_2))
        return -1e9

    return localSimVal + globalSimVal

simFunc = c.SIM_FUNC(CombinedSimFunc)

pAmino_1, seqStr_1, len_1 = c.PDB2CAminoPtrArray(len(aminoNVtr_1))
pAmino_2, seqStr_2, len_2 = c.PDB2CAminoPtrArray(len(aminoNVtr_2))
revTracePath = c.NW.Needleman_Wunsch(pAmino_1, pAmino_2, simFunc, 50, 0, len_seq_1=len_1, len_seq_2=len_2)

resStr_1, resStr_2 = c.NW.TracePath2AlignedStr(revTracePath, seqStr_1.decode('ascii'), seqStr_2.decode('ascii'))
res_1, res_2 = c.NW.TracePath2OnlyAlignedList(revTracePath, resis_1, resis_2)

rotMat, transVtr, a, b = scr.TMscoreAligned(res_1, res_2, len_2)

theChain = ph.GetChain(pdbFile1, 'A')
transAtomList = []
for residue in theChain.residues:
    for atom in residue.atoms.values():
        transAtomList.append(atom.transform(rotMat.T, transVtr).infoList)

ph.SavePDBFile('pdbs/1MBA_rotated.pdb', [ph.GetChain(transAtomList, 'A')])
#ph.SavePDBFile('pdbs/1atz_rotated.pdb', [ph.GetChain(transAtomList, 'A')])

print()
print(resStr_1)
print(resStr_2)
