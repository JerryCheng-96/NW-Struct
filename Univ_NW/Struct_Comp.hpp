using namespace std;

/*
 *  Struct Definitions
 */

typedef struct vector {
    double x;
    double y;
    double z;
} Vtr;

typedef struct amino_acid {
    Vtr* vector;
    char ss;
} Amino;

/*
 *  Function declarations
 */

extern "C" {
    double Vtr_length(Vtr* vtr);
    double Vtr_angle_cosine(Vtr* vtr_1, Vtr* vtr_2);
    double Sim_amino_vtrCos(void* elem_1, void* elem_2);

    int Show_Amino(Amino*);
}


double Vtr_length(Vtr* vtr) {
    //printf("x = %5f, y = %5f, z = %5f\n", vtr->x, vtr->y, vtr->z);
    return sqrt(vtr->x * vtr->x + vtr->y * vtr->y + vtr->z * vtr->z);
}

double Vtr_angle_cosine(Vtr* vtr_1, Vtr* vtr_2) {
    return (vtr_1->x * vtr_2->x + vtr_1->y * vtr_2->y + vtr_1->z * vtr_2->z) / (Vtr_length(vtr_1) * Vtr_length(vtr_2));
}

double Sim_amino_vtrCos(void* elem_1, void* elem_2) {
    Amino* amino_1 = (Amino*) elem_1;
    Amino* amino_2 = (Amino*) elem_2;

    // If the Sec Struct is different, assign a very negative similarity value
//    if (amino_1->ss == amino_2->ss)
//        return 1e2;

    // Calculating the cosine value of the angle
    return acos(Vtr_angle_cosine(amino_1->vector, amino_2->vector)) * (180 / acos(-1.0));
}

