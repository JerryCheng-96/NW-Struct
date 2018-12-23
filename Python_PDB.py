import math
import numpy as np
from Bio.PDB.PDBParser import PDBParser


## The "infinitesimal"
epsilon = 1e-9


def Len_Vtr(vtr):
    return math.sqrt(sum([comp * comp for comp in list(vtr)]))


def Calc_diheral(atoms_vtr):
    if not len(atoms_vtr) == 4:
        print("In Calc_diheral: wrong number of atoms!! Requires 4 atoms.")
        raise ValueError()
    vtr_10 = atoms_vtr[1] - atoms_vtr[0]
    vtr_12 = atoms_vtr[1] - atoms_vtr[2]
    vtr_23 = atoms_vtr[2] - atoms_vtr[3]
    
    norm_p1 = np.cross(list(vtr_10), list(vtr_12))
    norm_p2 = np.cross(list(-vtr_12), list(vtr_23))

    if np.dot(list(vtr_23), norm_p1) > 0:
        return math.degrees(math.acos(np.dot(norm_p1, norm_p2) / (Len_Vtr(norm_p1) * Len_Vtr(norm_p2))))
    else:
        return -math.degrees(math.acos(np.dot(norm_p1, norm_p2) / (Len_Vtr(norm_p1) * Len_Vtr(norm_p2))))


def Calc_phipsi(resis):
    phipsi_list = []
    
    ## At the very first resi...
    phipsi_list.append([0, Calc_diheral([resis[0]["N"].get_vector(), resis[0]["CA"].get_vector(), resis[0]["C"].get_vector(), resis[1]["N"].get_vector()])])
    
    ## Iterating the residues in between...
    for i in range(1, len(resis) - 1):
        phipsi_list.append([Calc_diheral([resis[i-1]["C"].get_vector(), resis[i]["N"].get_vector(), resis[i]["CA"].get_vector(), resis[i]["C"].get_vector()]),\
            Calc_diheral([resis[i]["N"].get_vector(), resis[i]["CA"].get_vector(), resis[i]["C"].get_vector(), resis[i+1]["N"].get_vector()])])

    ## The last residue...
    i = len(resis) - 1
    phipsi_list.append([Calc_diheral([resis[i-1]["C"].get_vector(), resis[i]["N"].get_vector(), resis[i]["CA"].get_vector(), resis[i]["C"].get_vector()]), 0])

    return phipsi_list


def Rodrigues_rot(vtr, axis, angle):
    return math.cos(angle) * vtr + np.cross(axis, vtr) * math.sin(angle) + \
        np.dot(axis, vtr) * (1 - math.cos(angle)) * axis


def SolveFor_rot(from_vtr, to_vtr):
    if not Len_Vtr(from_vtr) - Len_Vtr(to_vtr) <= epsilon:
        print("The two vectors should be of the same length.")
        raise ValueError()
    axis = np.cross(from_vtr, to_vtr)
    axis = axis / Len_Vtr(axis)     ## ... so the length of the vector is 1
    cos_angle = np.dot(from_vtr, to_vtr) / (Len_Vtr(from_vtr) * Len_Vtr(to_vtr))
    
    return (axis, math.acos(cos_angle))


def Get_AnimoNVtr(resis, window = 1):
    listNVtr = [np.array(list(resis[1]["N"].get_vector() - resis[0]["N"].get_vector()))]

    for i in range(1, len(resis) - 1):
        vtr_NC = list(resis[i]["N"].get_vector() - resis[i-1]["C"].get_vector())
        vtr_NCa = list(resis[i]["N"].get_vector() - resis[i]["CA"].get_vector())
        norm_vtr = np.cross(vtr_NC, vtr_NCa)
        norm_vtr = norm_vtr / Len_Vtr(norm_vtr)
        axis, angle = SolveFor_rot(norm_vtr, np.array([0,0,1]))
        oriNVtr = resis[i]["N"].get_vector() - resis[i+1]["N"].get_vector()
        listNVtr.append(Rodrigues_rot(np.array(list(oriNVtr)), axis, angle))
    
    listNVtr.append(np.array([0,0,0]))
    return listNVtr


if __name__ == "__main__":
    pp = PDBParser(PERMISSIVE=1)

    id = "1MBA"
    filename = "C:/Users/Jerry/CLionProjects/PDB2RMSD_2/example/1MBA.pdb"
    struct = pp.get_structure(id, filename)

    model = struct[0]
    chain = model["A"]

    ## Getting the animo acids residues
    resis = [resi for resi in chain if resi.get_id()[0] == ' ']
    phipsi_list = Calc_phipsi(resis)

    i = 0
    for line in phipsi_list:
        print(str(i) + " " + str(line))
        i += 1

