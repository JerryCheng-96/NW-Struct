from ctypes import *
import Python_PDB as ppdb
libstruct = CDLL("Univ_NW/libstruct.so")


# Structures in Struct_Comp.hpp
class CVtr(Structure):
    _fields_ = [("x", c_double), ("y", c_double), ("z", c_double)]


class CAmino(Structure):
    _fields_ = [("vector", POINTER(CVtr)), ("ss", c_char)]


# Functions in Struct_Comp.hpp
def CVtrLength(CVtr):       # Returning the length (norm) of a CVtr
    libstruct.Vtr_length.restype = c_double
    return libstruct.Vtr_length(CVtr)


def CVtrAngleCosine(CVtr_1, CVtr_2):    # Returning the cosine value of the angle between two vectors
    libstruct.Vtr_angle_cosine.restype = c_double
    return libstruct.Vtr_angle_cosine(CVtr_1, CVtr_2)


# Python helper functions
def Np2CVtr(np_arr):
    return CVtr(np_arr[0], np_arr[1], np_arr[2])


def PDB2CAminoPtrArray(filename):
    nVtrList = ppdb.Get_AnimoNVtr(ppdb.Get_AminoResidues(filename))
    dsspList = ppdb.Get_DSSP(filename)
    assert len(nVtrList) == len(dsspList)
    seq_len = len(nVtrList)

    pAminoArray = (POINTER(CAmino) * seq_len)()
    pAminoArray = cast(pAminoArray, POINTER(POINTER(CAmino)))

    aminos_seq = ''

    for i in range(0, seq_len):
        pAminoArray[i] = pointer(CAmino(pointer(Np2CVtr(nVtrList[i])), ord(dsspList[i][2])))
        aminos_seq += dsspList[i][1]

    return pAminoArray, bytes(aminos_seq, encoding="ASCII"), seq_len

