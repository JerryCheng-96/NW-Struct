import re
from Vector_Matrix import *
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBParser import PDBParser


def Get_PDBStruct(filename):
    pp = PDBParser(PERMISSIVE=1)
    try:
        id = re.findall(r"([^/]+?)\.pdb$", filename)[0]
    except IndexError as ie:
        print("The input file should be a .pdb.")
        raise(ie)

    return pp.get_structure(id, filename)[0], filename


def Get_DSSP(filename, struct_no=0):
    pp = PDBParser(PERMISSIVE=1)
    try:
        id = re.findall(r"([^/]+?)\.pdb$", filename)[0]
    except IndexError as ie:
        print("The input file should be a .pdb.")
        raise (ie)

    struct = pp.get_structure(id, filename)[0], filename

    try:
        model = struct[struct_no]
    except IndexError as ie:
        print("No such structure with number \"" + struct_no + "\"")
        raise(ie)

    dssp = DSSP(model, filename, dssp="mkdssp")
    return [dssp[i] for i in list(dssp.keys())]


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


def Get_AminoResidues(struct, struct_no=0, chain_id="A"):
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
