from Vector_Matrix import *


class Atom:
    def __init__(self, infoList):
        self.infoList = infoList

    @property
    def vector(self):
        return np.array([float(val) for val in self.infoList[6:9]])

    @property
    def type(self):
        return self.infoList[2]

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

    @property
    def AtomsList(self):
        atomsList = []
        for residue in self.residues:
            atomsList += residue.atoms.values()
        return atomsList

    @property
    def CAAtomsList(self):
        return [atom for atom in self.AtomsList if atom.type == 'CA']

    def __len__(self):
        return len(self.residues)


    # TODO: Implement the window function
    def Get_AminoNVtr(self, window = 1):
        listNVtr = [self.residues[1].atoms['N'].vector - self.residues[0].atoms['N'].vector]

        for i in range(1, len(self.residues) - 1):
            vtr_NC = self.residues[i].atoms['N'].vector - self.residues[i - 1].atoms['C'].vector
            vtr_NCa = self.residues[i].atoms['N'].vector - self.residues[i].atoms['CA'].vector
            norm_vtr = np.cross(vtr_NC, vtr_NCa)
            norm_vtr = norm_vtr / Len_Vtr(norm_vtr)
            axis, angle = SolveFor_rot(norm_vtr, np.array([0,0,1]))
            oriNVtr = self.residues[i].atoms['N'].vector - self.residues[i + 1].atoms['N'].vector
            listNVtr.append(Rodrigues_rot(oriNVtr, axis, angle))

        listNVtr.append(np.array([0,0,0]))
        return listNVtr

    def Get_CAContactMap(self):
        contactMap = np.empty((len(self), len(self)))
        for i in range(0, len(self)):
            for j in range(i, len(self)):
                contactMap[i][j] = Len_Vtr(self.residues[i].atoms['CA'].vector - self.residues[j].atoms['CA'].vector)
        return contactMap + contactMap.T


    def Get_SurroundVectorSet(self, low_cutoff, high_cutoff, keep=10):
        neighborAminos = []

        for i in range(0, len(self)):
            residueLenAminos = {}

            vtr_CaN = self.residues[i].atoms['CA'].vector - self.residues[i - 1].atoms['N'].vector
            vtr_CaC = self.residues[i].atoms['CA'].vector - self.residues[i].atoms['C'].vector
            norm_vtr = np.cross(vtr_CaN, vtr_CaC)
            norm_vtr = norm_vtr / Len_Vtr(norm_vtr)
            axis, angle = SolveFor_rot(norm_vtr, np.array([0,0,1]))

            for j in range(0, len(self)):
                theVtr = self.residues[i].atoms['CA'].vector - self.residues[j].atoms['CA'].vector
                theVtrLen = Len_Vtr(theVtr)
                if low_cutoff < theVtrLen < high_cutoff:
                    residueLenAminos[Len_Vtr(theVtr)] = theVtr

            residueSurroundVtr = [Rodrigues_rot(residueLenAminos[aaLen], axis, angle)\
                                  for aaLen in sorted(residueLenAminos.keys())][:keep]
            neighborAminos.append(residueSurroundVtr)

        return neighborAminos


def Get_NeighborAminoNo(contactMap, low_cutoff, high_cutoff):
    neighborVector = []

    for i in range(0, len(contactMap)):
        theResidueNeighbor = []
        for j in range(0, len(contactMap)):
            if low_cutoff < contactMap[i][j] <= high_cutoff:
                theResidueNeighbor.append(j)
        neighborVector.append(theResidueNeighbor)

    return neighborVector


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
#    for i in range(1, 395):
#        try:
#            pdbAminoDict[i] == None
#        except KeyError:
#            print(i)

    return Chain([pdbAminoDict[aminoNo] for aminoNo in sorted(pdbAminoDict.keys())])


if __name__ == "__main__":
    pdbAtomsList = ReadPDBAsAtomsList("pdbs/Legacy/4xt3.pdb")
    chainA = GetChain(pdbAtomsList, 'A')
    chainA.Get_AminoNVtr()
    ##cm = chainA.Get_CAContactMap()
    ##nv = Get_NeighborAminoNo(cm, 4.0, 12.0)
    neiVtr = chainA.Get_SurroundVectorSet(3.8, 12, 10)
    simValArray = []
    for i in range(90, 110):
        simValArray.append(Calc_SimVectorSet(neiVtr[100], neiVtr[i]))
    SavePDBFile("pdbs/4xt3_new.pdb", [chainA])
    print()
