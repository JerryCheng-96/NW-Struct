#include <iostream>
#include <fstream>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "Univ_NW_Debug.hpp"
#include "Univ_NW.hpp"
#include "Blosum_Helper.hpp"
#include "Struct_Comp.hpp"

using namespace std;

extern "C" {
    float relay_func(void*, void*);
    int funcptr_test();
    int test(int (*funcptr)());
    int test_array(Vtr* [], int);
    int test_str(char*);
    int Show_Amino(Amino*);

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


float VoidPtrTest(void* aa_1, void* aa_2) {
    ;
}

