from ctypes import *
import PDB_Parse as pdbPar
libc = CDLL("Univ_NW/libstruct.so")


# ****************
# ***VtrMat.hpp***
# ****************

# Structs
class CVtr(Structure):
    _fields_ = [("x", c_double), ("y", c_double), ("z", c_double)]


# Functions
def CVtrLength(CVtr):       # Returning the length (norm) of a CVtr
    libc.Vtr_length.restype = c_double
    return libc.Vtr_length(CVtr)


def CVtrAngleCosine(CVtr_1, CVtr_2):    # Returning the cosine value of the angle between two vectors
    libc.Vtr_angle_cosine.restype = c_double
    return libc.Vtr_angle_cosine(CVtr_1, CVtr_2)


# Python helper functions
def Np2CVtr(np_arr):
    return CVtr(np_arr[0], np_arr[1], np_arr[2])


# *********************
# ***Struct_Comp.cpp***
# *********************

# Structs
class CAmino(Structure):
    _fields_ = [("vector", POINTER(CVtr)), ("ss", c_char)]


# Python helper functions
##def PDB2CAminoPtrArray(nVtrList, dsspList):
##    assert (nVtrList is not None and dsspList is not None)
##    assert len(nVtrList) == len(dsspList)
##    seq_len = len(nVtrList)
##    #nVtrList = pdbPar.Get_AnimoNVtr(pdbPar.Get_AminoResidues(struct))
##    #dsspList = pdbPar.Get_DSSP(struct, filename)
##
##    pAminoArray = (POINTER(CAmino) * seq_len)()
##    pAminoArray = cast(pAminoArray, POINTER(POINTER(CAmino)))
##
##    aminos_seq = ''
##
##    for i in range(0, seq_len):
##        pAminoArray[i] = pointer(CAmino(pointer(Np2CVtr(nVtrList[i])), ord(dsspList[i][2])))
##        aminos_seq += dsspList[i][1]
##
##    return pAminoArray, bytes(aminos_seq, encoding="ASCII"), seq_len


def PDB2CAminoPtrArray(nVtrList):
    assert (nVtrList is not None)
    seq_len = len(nVtrList)

    pAminoArray = (POINTER(CAmino) * seq_len)()
    pAminoArray = cast(pAminoArray, POINTER(POINTER(CAmino)))

    aminos_seq = ''

    for i in range(0, seq_len):
        pAminoArray[i] = pointer(CAmino(pointer(Np2CVtr(nVtrList[i])), ord('-')))
        aminos_seq += 'A'

    return pAminoArray, bytes(aminos_seq, encoding="ASCII"), seq_len



# *****************
# ***Univ_NW.hpp***
# *****************

# Defining a class called NW
class NW:
    @staticmethod
    def Needleman_Wunsch(seq_1, seq_2, sim_func, gap, ext, len_seq_1=None, len_seq_2=None):
        if len_seq_1 is None:
            len_seq_1 = len(seq_1)
        if len_seq_2 is None:
            len_seq_2 = len(seq_2)

        libc.NW_Align.restype = c_char_p;
        tracePath = libc.NW_Align(seq_1, len_seq_1, seq_2, len_seq_2, sim_func, gap, ext)
        return tracePath.decode('ascii')

    @staticmethod
    def TracePath2AlignedList(revTracePath, seq_1, seq_2):
        tracePath = revTracePath[::-1]
        res_1 = []
        res_2 = []
        idx_1 = 0
        idx_2 = 0

        for i in tracePath:
            if i == '.':
                res_1.append(seq_1[idx_1])
                res_2.append(seq_2[idx_2])
                idx_1 += 1
                idx_2 += 1
            elif i == '-':
                res_1.append(seq_1[idx_1])
                res_2.append(None)
                idx_1 += 1
            elif i == '+':
                res_2.append(seq_2[idx_2])
                res_1.append(None)
                idx_2 += 1

        return res_1, res_2

    @staticmethod
    def TracePath2OnlyAlignedList(revTracePath, seq_1, seq_2):
        tracePath = revTracePath[::-1]
        res_1 = []
        res_2 = []
        idx_1 = 0
        idx_2 = 0

        for i in tracePath:
            if i == '.':
                res_1.append(seq_1[idx_1])
                res_2.append(seq_2[idx_2])
                idx_1 += 1
                idx_2 += 1
            elif i == '-':
                idx_1 += 1
            elif i == '+':
                idx_2 += 1

        return res_1, res_2

    @staticmethod
    def TracePath2AlignedStr(revTracePath, seqStr_1, seqStr_2):
        tracePath = revTracePath[::-1]
        resStr_1 = ""
        resStr_2 = ""
        idx_1 = 0
        idx_2 = 0

        for i in tracePath:
            if i == '.':
                resStr_1 += seqStr_1[idx_1]
                resStr_2 += seqStr_2[idx_2]
                idx_1 += 1
                idx_2 += 1
            elif i == '-':
                resStr_1 += seqStr_1[idx_1]
                resStr_2 += '-'
                idx_1 += 1
            elif i == '+':
                resStr_2 += seqStr_2[idx_2]
                resStr_1 += '-'
                idx_2 += 1

        return resStr_1, resStr_2


# Similarity functions
def SimAminoCosine(ptr_1, ptr_2):
    ptr_1 = cast(ptr_1, POINTER(CAmino))
    ptr_2 = cast(ptr_2, POINTER(CAmino))
    return 100 * CVtrAngleCosine(ptr_1.contents.vector, ptr_2.contents.vector)

SIM_FUNC = CFUNCTYPE(c_float, c_void_p, c_void_p)
cos_simFunc = SIM_FUNC(SimAminoCosine)

