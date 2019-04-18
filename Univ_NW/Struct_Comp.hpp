using namespace std;

/*
 *  Struct Definitions
 */

typedef struct amino_acid {
    Vtr* vector;
    char ss;
} Amino;


typedef struct pyAminoAcid {
    int num;
} PyAmino;


/*
 *  Function declarations
 */

extern "C" {
    double Sim_amino_vtrCos(void* elem_1, void* elem_2);
    int Show_Amino(Amino*);
}


/*
 *  Function definitions
 */

double Sim_amino_vtrCos(void* elem_1, void* elem_2) {
    Amino* amino_1 = (Amino*) elem_1;
    Amino* amino_2 = (Amino*) elem_2;

    // If the Sec Struct is different, assign a very negative similarity value
//    if (amino_1->ss == amino_2->ss)
//        return 1e2;

    // Calculating the cosine value of the angle
//    return acos(Vtr_angle_cosine(amino_1->vector, amino_2->vector)) * (180 / acos(-1.0));
    return 10 * Vtr_angle_cosine(amino_1->vector, amino_2->vector);
}

