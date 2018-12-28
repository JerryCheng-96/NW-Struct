import re
from ctypes import *
libstruct = CDLL("Univ_NW/libstruct.so")


gBlosumMat = {}


def Read_Blosum(filename):
    f = open(filename)
    blosumLines = f.readlines()
    firstRow = [' ']
    firstRow += (re.split(r"\s*", blosumLines[0].strip()))
    blosumDict = {}

    for line in blosumLines[1:]:
        line = re.split(r"\s*", line.strip())
        blosumDict[line[0]] = {}
        for i in range(1, len(line)):
            blosumDict[line[0]][firstRow[i]] = int(line[i])

    return blosumDict


def BlosumSimFunc(aa_1, aa_2):
    aa_1 = cast(aa_1, c_char_p)
    aa_2 = cast(aa_2, c_char_p)
    aa_1 = aa_1.value.decode("ASCII")[0]
    aa_2 = aa_2.value.decode("ASCII")[0]

    return gBlosumMat[aa_1][aa_2]


def TracePath2Align(revTracePath, seqStr_1, seqStr_2):
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

    return (resStr_1, resStr_2)


gBlosumMat = Read_Blosum("Univ_NW/BLOSUM62.txt")
if __name__ == "__main__":
    SIM_FUNC = CFUNCTYPE(c_float, c_void_p, c_void_p)
    blosum_sim_func = SIM_FUNC(BlosumSimFunc)

#    seq_1 = b'SLSAAEADLAGKSWAPVFANKNANGLDFLVALFEKFPDSANFFADFKGKSVADIKASPKLRDVSSRIFTRLNEFVNNAANAGKMSAMLSQFAKEHVGFGVGSAQFENVRSMFPGFVASVAAPPAGADAAWTKLFGLIIDALKAAGA'
#    seq_2 = b'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG'
    seq_1 = b'MAAAAASLRRTVLGPRGVGLPGASAPGLLGGARSRQLPLRTPQAVSLSSKSGPSRGRKVMLSALGMLAAGGAGLAVALHSAVSASDLELHPPSYPWSHRGLLSSLDHTSIRRGFQVYKQVCSSCHSMDYVAYRHLVGVCYTEEEAKALAEEVEVQDGPNDDGEMFMRPGKLSDYFPKPYPNPEAARAANNGALPPDLSYIVRARHGGEDYVFSLLTGYCEPPTGVSLREGLYFNPYFPGQAIGMAPPIYTEVLEYDDGTPATMSQVAKDVATFLRWASEPEHDHRKRMGLKMLLMMGLLLPLTYAMKRHKWSVLKSRKLAYRPPK'
    seq_2 = b'MLSVAARSGPFAPVLSATSRGVAGALRPLVQAAVPATSESPVLDLKRSVLCRESLRGQAAGRPLVASVSLNVPASVRYSHTDIKVPDFSDYRRPEVLDSTKSSKESSEARKGFSYLVTATTTVGVAYAAKNVVSQFVSSMSASADVLAMSKIEIKLSDIPEGKNMAFKWRGKPLFVRHRTKKEIDQEAAVEVSQLRDPQHDLERVKKPEWVILIGVCTHLGCVPIANAGDFGGYYCPCHGSHYDASGRIRKGPAPLNLEVPSYEFTSDDMVIVG'

    libstruct.Blosum_NW_Relay.restype = c_char_p
    tracePath = libstruct.Blosum_NW_Relay(seq_1, len(seq_1), seq_2, len(seq_2), blosum_sim_func, 11, 1)
    res_1, res_2 = TracePath2Align(tracePath.decode('ascii'), seq_1.decode('ascii'), seq_2.decode('ascii'))
    print(res_1)
    print(tracePath.decode('ascii')[::-1])
    print(res_2)
