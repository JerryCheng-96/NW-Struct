import re
import numpy as np
from Vector_Matrix import *


class Atom:
    def __init__(self, infoList):
        self.infoList = infoList

    @property
    def vector(self):
        return np.array([float(val) for val in self.infoList[6:9]])

    @property
    def n(self):
        return int(self.infoList[1])

    def transform(self, rotMat, transVtr):
        vtr = np.dot(self.vector, rotMat) + transVtr
        self.infoList[6] = '%.3f' % vtr[0]
        self.infoList[7] = '%.3f' % vtr[1]
        self.infoList[8] = '%.3f' % vtr[2]
        return self


class Residue:
    def __init__(self, atomDict):
        self.atoms = atomDict


class Chain:
    def __init__(self, residueList):
        self.residues = residueList


def ReadPDBAsAtomsList(filename):
    f = open(filename, 'r')
    content = f.readlines()
    atomLines = [theLine for theLine in content if theLine[:4] == "ATOM"]
    pdbAtomsList = []

    for line in atomLines:
        theAtom = []
        theAtom.append(line[:4])
        theAtom.append(line[4:11].strip())
        theAtom.append(line[11:17].strip())
        theAtom.append(line[17:20].strip())
        theAtom.append(line[20:22].strip())
        theAtom.append(line[22:26].strip())
        theAtom.append(line[26:38].strip())
        theAtom.append(line[38:46].strip())
        theAtom.append(line[46:54].strip())
        theAtom.append(line[54:60].strip())
        theAtom.append(line[60:66].strip())
        theAtom.append(line[66:80].strip())
        pdbAtomsList.append(theAtom)

    return pdbAtomsList


def FormatAtomLine(atomInfoList):
    return '%s%7s'%(atomInfoList[0], atomInfoList[1]) + 2 * ' ' + atomInfoList[2].ljust(4) + atomInfoList[3].ljust(4) + \
           atomInfoList[4] + atomInfoList[5].rjust(4) + atomInfoList[6].rjust(12) + atomInfoList[7].rjust(8) + \
           atomInfoList[8].rjust(8) + atomInfoList[9].rjust(6) + atomInfoList[10].rjust(6) + ' ' * 11 + \
           atomInfoList[11].ljust(3) + '\n'

def SavePDBFile(filename, chainList):
    f = open(filename, 'w')
    pdbFileLines = []
    atomsDict = {}

    for chain in chainList:
        for residue in chain.residues:
            for atom in residue.atoms.values():
                atomsDict[atom.n] = atom
        pdbFileLines += [FormatAtomLine(atomsDict[atomNo].infoList) + '\n' for atomNo in sorted(atomsDict.keys())]
        pdbFileLines += "TER\n"

    pdbFileLines += 'END\n'
    f.writelines(pdbFileLines)
    f.flush()
    f.close()


def GetChain(pdbAtomsList, chainId):
    atomsList = [atom for atom in pdbAtomsList if atom[4] == chainId]
    pdbAminoDict = {}

    for atom in atomsList:
        try:
            pdbAminoDict[int(atom[5])].atoms[atom[2]] = Atom(atom)
        except KeyError:
            pdbAminoDict[int(atom[5])] = Residue({})
            pdbAminoDict[int(atom[5])].atoms[atom[2]] = Atom(atom)
##    for i in range(1, 395):
##        try:
##            pdbAminoDict[i] == None
##        except KeyError:
##            print(i)

    return Chain([pdbAminoDict[aminoNo] for aminoNo in sorted(pdbAminoDict.keys())])


def Get_AminoNVtr(chain, window = 1):
    listNVtr = [chain.residues[1].atoms['N'].vector - chain.residues[0].atoms['N'].vector]

    for i in range(1, len(chain.residues) - 1):
        vtr_NC = chain.residues[i].atoms['N'].vector - chain.residues[i-1].atoms['C'].vector
        vtr_NCa = chain.residues[i].atoms['N'].vector - chain.residues[i].atoms['CA'].vector
        norm_vtr = np.cross(vtr_NC, vtr_NCa)
        norm_vtr = norm_vtr / Len_Vtr(norm_vtr)
        axis, angle = SolveFor_rot(norm_vtr, np.array([0,0,1]))
        oriNVtr = chain.residues[i].atoms['N'].vector - chain.residues[i+1].atoms['N'].vector
        listNVtr.append(Rodrigues_rot(oriNVtr, axis, angle))

    listNVtr.append(np.array([0,0,0]))
    return listNVtr


if __name__ == "__main__":
    pdbAtomsList = ReadPDBAsAtomsList("pdbs/Legacy/4xt3.pdb")
    chainA = GetChain(pdbAtomsList, 'A')
    Get_AminoNVtr(chainA)
    SavePDBFile("pdbs/4xt3_new.pdb", [chainA])
    print()



