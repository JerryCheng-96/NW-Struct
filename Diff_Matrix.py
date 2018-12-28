import re

def ReadInMatrix(filename):
    f = open(filename)
    theLines = f.readlines()
    theLinesLists = []
    for line in theLines:
        theLinesLists.append(re.split("\s*", line))
    return theLinesLists

def DiffMatrix(mat_1, mat_2):
    for i in range(0, min(len(mat_1), len(mat_2))):
        for j in range(0, min(len(mat_1[i]), len(mat_2[i]))):
            if mat_1[i][j] != mat_2[i][j]:
                print("Diff at (" + str(i) + ", " + str(j) + "): Mat_1 = " + str(mat_1[i][j]) + ", Mat_2 = " + str(mat_2[i][j]))

DiffMatrix(ReadInMatrix("output_cao.txt"), ReadInMatrix("output_trace.txt"))
