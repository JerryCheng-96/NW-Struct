#include <iostream>
#include <fstream>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "Univ_NW.hpp"
#include "Struct_Comp.hpp"

using namespace std;

extern "C" {
    float relay_func(void*, void*);
    int funcptr_test();
    int test(int (*funcptr)());
    int test_array(Vtr* [], int);
    int test_str(char*);
    int Show_Amino(Amino*);

    int Blosum_NW_Relay(char*, int, char*, int, float(*sim_func)(void*, void*), int, int);
    float VoidPtrTest(void*, void*);
}

int main() {

    return 0;
}

float relay_func(void* elem_1, void* elem_2) {
    return 0;
}

int funcptr_test() {
    printf("Now we are in funcptr_test!!!\n");
}

int test(int (*funcptr)()) {
    funcptr();
    printf("Hello world!\n");
    return 0;
}

int test_array(Vtr* vtrArray[], int lenArray) {
    for (int i = 0; i < lenArray; i++) {
        Vtr* vtr = vtrArray[i];
        printf("Now at Vtr[%d], x = %5f, y = %5f, z = %5f\n", i, vtr->x, vtr->y, vtr->z);
    }
    return 0;
}

int test_str(char* str) {
    printf("%s\n", str);
    return 0;
}

int Show_Amino(Amino* ami) {
    Vtr* vtr = ami->vector;
    printf("Vtr = (%5f, %5f, %5f)\t", vtr->x, vtr->y, vtr->z);
    printf("%c\n", ami->ss);

    return 0;
}

int Blosum_NW_Relay(char* chr_seq_1, int len_seq_1, char* chr_seq_2, int len_seq_2,
                    float (*sim_func)(void*, void*), int gap_open, int gap_ext) {
    void** seq_1 = (void**)malloc(sizeof(void*) * (len_seq_1 + 10));
    void** seq_2 = (void**)malloc(sizeof(void*) * (len_seq_2 + 10));

    for (int i = 0; i < len_seq_1; i++) {
        seq_1[i] = (void*)(chr_seq_1 + i);
    }

    for (int i = 0; i < len_seq_2; i++) {
        seq_2[i] = (void*)(chr_seq_2 + i);
    }

    NW_Align(seq_1, chr_seq_1, len_seq_1, seq_2, chr_seq_2, len_seq_2, sim_func, gap_open, gap_ext);
}

float VoidPtrTest(void* aa_1, void* aa_2) {
    ;
}

