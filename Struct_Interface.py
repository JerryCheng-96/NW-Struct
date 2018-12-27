from ctypes import *
libstruct = CDLL("/mnt/c/Users/Jerry/Desktop/Grad Design/Univ_NW/libstruct.so")


# Corresponding structures
class CVtr(Structure):
    _fields_ = [("x", c_double), ("y", c_double), ("z", c_double)]


class CAmino(Structure):
    _fields_ = [("vector", POINTER(CVtr)), ("ss", c_char)]


# Functions
def CVtrLength(CVtr):
    libstruct.Vtr_length.restype = c_double
    return libstruct.Vtr_length(CVtr)


def CVtrAngleCosine(CVtr_1, CVtr_2):
    libstruct.Vtr_angle_cosine.restype = c_double
    return libstruct.Vtr_angle_cosine(CVtr_1, CVtr_2)
