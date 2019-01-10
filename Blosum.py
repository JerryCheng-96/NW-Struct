import re
import C_Interface as C

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
    aa_1 = C.cast(aa_1, C.c_char_p)
    aa_2 = C.cast(aa_2, C.c_char_p)
    aa_1 = aa_1.value.decode("ASCII")[0]
    aa_2 = aa_2.value.decode("ASCII")[0]

    return gBlosumMat[aa_1][aa_2]


gBlosumMat = Read_Blosum("Univ_NW/BLOSUM62.txt")
if __name__ == "__main__":
    SIM_FUNC = C.CFUNCTYPE(C.c_float, C.c_void_p, C.c_char_p)
    blosum_sim_func = SIM_FUNC(BlosumSimFunc)

#    seq_1 = b'SLSAAEADLAGKSWAPVFANKNANGLDFLVALFEKFPDSANFFADFKGKSVADIKASPKLRDVSSRIFTRLNEFVNNAANAGKMSAMLSQFAKEHVGFGVGSAQFENVRSMFPGFVASVAAPPAGADAAWTKLFGLIIDALKAAGA'
#    seq_2 = b'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG'
    seq_1 = b'MAAAAASLRRTVLGPRGVGLPGASAPGLLGGARSRQLPLRTPQAVSLSSKSGPSRGRKVMLSALGMLAAGGAGLAVALHSAVSASDLELHPPSYPWSHRGLLSSLDHTSIRRGFQVYKQVCSSCHSMDYVAYRHLVGVCYTEEEAKALAEEVEVQDGPNDDGEMFMRPGKLSDYFPKPYPNPEAARAANNGALPPDLSYIVRARHGGEDYVFSLLTGYCEPPTGVSLREGLYFNPYFPGQAIGMAPPIYTEVLEYDDGTPATMSQVAKDVATFLRWASEPEHDHRKRMGLKMLLMMGLLLPLTYAMKRHKWSVLKSRKLAYRPPK'
    seq_2 = b'MLSVAARSGPFAPVLSATSRGVAGALRPLVQAAVPATSESPVLDLKRSVLCRESLRGQAAGRPLVASVSLNVPASVRYSHTDIKVPDFSDYRRPEVLDSTKSSKESSEARKGFSYLVTATTTVGVAYAAKNVVSQFVSSMSASADVLAMSKIEIKLSDIPEGKNMAFKWRGKPLFVRHRTKKEIDQEAAVEVSQLRDPQHDLERVKKPEWVILIGVCTHLGCVPIANAGDFGGYYCPCHGSHYDASGRIRKGPAPLNLEVPSYEFTSDDMVIVG'

    C.libc.Blosum_NW_Relay.restype = C.c_char_p
    tracePath = C.libc.Blosum_NW_Relay(seq_1, len(seq_1), seq_2, len(seq_2), blosum_sim_func, 11, 1)
    res_1, res_2 = C.NW.TracePath2AlignedStr(tracePath.decode('ascii'), seq_1.decode('ascii'), seq_2.decode('ascii'))
    print(res_1)
    print(tracePath.decode('ascii')[::-1])
    print(res_2)
