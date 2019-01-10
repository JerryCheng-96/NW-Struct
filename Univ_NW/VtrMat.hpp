
typedef struct vector {
    double x;
    double y;
    double z;
} Vtr;

extern "C" {
    double Vtr_length(Vtr* vtr);
    double Vtr_angle_cosine(Vtr* vtr_1, Vtr* vtr_2);
}

double Vtr_length(Vtr* vtr) {
    return sqrt(vtr->x * vtr->x + vtr->y * vtr->y + vtr->z * vtr->z);
}

double Vtr_angle_cosine(Vtr* vtr_1, Vtr* vtr_2) {
    return (vtr_1->x * vtr_2->x + vtr_1->y * vtr_2->y + vtr_1->z * vtr_2->z) / (Vtr_length(vtr_1) * Vtr_length(vtr_2));
}