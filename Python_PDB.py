import re
import math
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP


# The "infinitesimal"
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


def Kabsch(from_points, to_points):
    ## Checking input...
    assert len(from_points) == len(to_points)

    ## Centering the points sets
    centroid_from = np.sum(from_points, axis=0) / len(from_points)
    centroid_to = np.sum(to_points, axis=0) / len(to_points)
    from_points = from_points.astype(dtype=np.float64) - centroid_from
    to_points = to_points.astype(dtype=np.float64) - centroid_to

    ## Cross variance matrix cv_mat
    cv_mat = np.dot(from_points.T,  to_points)

    ## SVD
    mat_U, mat_S, mat_VT = np.linalg.svd(cv_mat)

    ## Solving for final rotation matrix mat_rot
    val_d = np.linalg.det(np.dot(mat_VT.T, mat_U.T))
    mat_rot = np.dot(mat_VT.T, np.array([[1, 0, 0],
                                         [0, 1, 0],
                                         [0, 0, val_d]]))
    mat_rot = np.dot(mat_rot, mat_U.T)

    return mat_rot


def Get_RotMatrix(roll=None, pitch=None, yaw=None):
    if roll is None:
        rot_xangle = np.random.rand()
    else:
        rot_xangle = roll
    if pitch is None:
        rot_yangle = np.random.rand()
    else:
        rot_yangle = pitch
    if yaw is None:
        rot_zangle = np.random.rand()
    else:
        rot_zangle = yaw

    rot_xmat = np.array([[1, 0, 0],
                         [0, math.cos(rot_xangle), -math.sin(rot_xangle)],
                         [0, math.sin(rot_xangle), math.cos(rot_xangle)]])
    rot_ymat = np.array([[math.cos(rot_yangle), 0, math.sin(rot_yangle)],
                         [0, 1, 0],
                         [-math.sin(rot_yangle), 0, math.cos(rot_yangle)]])
    rot_zmat = np.array([[math.cos(rot_zangle), -math.sin(rot_zangle), 0],
                         [math.sin(rot_zangle), math.cos(rot_zangle), 0],
                         [0, 0, 1]])

    rot_allmat = np.dot(np.dot(rot_zmat, rot_ymat), rot_xmat)

    return rot_allmat


def Get_DSSP(filename, struct_no=0):
    pp = PDBParser(PERMISSIVE=1)
    try:
        pdb_id = re.findall(r"([^/]+?)\.pdb$", filename)[0]
    except IndexError as ie:
        print("The input file should be a .pdb.")
        raise(ie)

    struct = pp.get_structure(pdb_id, filename)

    try:
        model = struct[struct_no]
    except IndexError as ie:
        print("No such structure with number \"" + struct_no + "\"")
        raise(ie)

    dssp = DSSP(model, filename, dssp="mkdssp")
    return [dssp[i] for i in list(dssp.keys())]


def Get_AminoResidues(filename, struct_no=0, chain_id="A"):
    pp = PDBParser(PERMISSIVE=1)
    try:
        id = re.findall(r"([^/]+?)\.pdb$", filename)[0]
    except IndexError as ie:
        print("The input file should be a .pdb.")
        raise(ie)

    struct = pp.get_structure(id, filename)

    model = struct[struct_no]
    chain = model[chain_id]

    # Getting the animo acids residues
    return [resi for resi in chain if resi.get_id()[0] == ' ']


if __name__ == "__main__":      # Test drive / codes
    rand_from_points = np.random.rand(15, 3)                                        # Generating a random set of 3-D points
    centroid_from = np.sum(rand_from_points, axis=0) / len(rand_from_points)
    rand_from_points = rand_from_points.astype(dtype=np.float64) - centroid_from    # Now centeroid at the coordinate origin
    rot_allmat = Get_RotMatrix()
    rand_to_points = [np.dot(rot_allmat, rand_from_points[0])]
    for i in range(1, len(rand_from_points)):
        rand_to_points.append(np.dot(rot_allmat, rand_from_points[i]))
    rand_to_points = np.array(rand_to_points)                                       # Rotating the points
    for i in range(0, 3):
        rand_to_points[np.random.randint(0, 14)] += \
            np.random.rand(1, 3)[0] * .1                                            # Adding pertubation

    kab_rotmat = Kabsch(rand_from_points, rand_to_points)                           # "Kabsching"...

    print("The randomly generated rotation matrix:")
    print(rot_allmat)
    print("\nThe Kabsch deduced rotation matrix:")
    print(kab_rotmat)
