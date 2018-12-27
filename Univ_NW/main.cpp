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