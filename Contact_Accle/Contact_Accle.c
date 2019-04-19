#include <math.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>


// Definition: Accessing the data at a certain position...
#define PTR_DOUBLE_ELEM_3(PDATA, STRD, I, J, K)     (double*)((PDATA) + (I) * (STRD)[0] + (J) * (STRD)[1] + (K) * STRD[2])
#define PTR_DOUBLE_ELEM_2(PDATA, STRD, I, J)        (double*)((PDATA) + (I) * (STRD)[0] + (J) * (STRD)[1])
#define PTR_DOUBLE_VTR(PVTR, STRD, I)               (double*)((PVTR) + (I) * (STRD))


// Definition: Simple vector operations
#define VTR_DOT_3D(PVTR1, STRD1, PVTR2, STRD2) (\
    (*(double*)((PVTR1) + (STRD1) * 0)) * (*(double*)((PVTR2) + (STRD2) * 0)) + \
    (*(double*)((PVTR1) + (STRD1) * 1)) * (*(double*)((PVTR2) + (STRD2) * 1)) + \
    (*(double*)((PVTR1) + (STRD1) * 2)) * (*(double*)((PVTR2) + (STRD2) * 2))   )
#define VTR_LEN_3D(PVTR, STRD) sqrt( \
    pow(*(double*)(PVTR), 2) + \
    pow(*(double*)((PVTR) + (STRD) * 1), 2) + \
    pow(*(double*)((PVTR) + (STRD) * 2), 2)   )
#define VTR_ANGLECOS(PVTR1, STRD1, PVTR2, STRD2) \
    VTR_DOT_3D((PVTR1), (STRD1), (PVTR2), (STRD2)) / ( \
        VTR_LEN_3D((PVTR1), (STRD1)) * \
        VTR_LEN_3D((PVTR2), (STRD2))   )

#define PRINT_VTR(PVTR, STRD) \
    printf("[%lf\t%lf\t%lf]\n", \
           *(double*)((PVTR) + 0 * (STRD)), \
           *(double*)((PVTR) + 1 * (STRD)), \
           *(double*)((PVTR) + 2 * (STRD))  )


/*
    Vector operations that require returning new vector
*/
PyObject* VTR_CROSS_3D(char* pVtr1, int strd1, char* pVtr2, int strd2) {
    double vtr1[] = {*(double*)((pVtr1) + (strd1) * 0), *(double*)((pVtr1) + (strd1) * 1), *(double*)((pVtr1) + (strd1) * 2)};
    double vtr2[] = {*(double*)((pVtr2) + (strd2) * 0), *(double*)((pVtr2) + (strd2) * 1), *(double*)((pVtr2) + (strd2) * 2)};

    //
    //for (int i = 0; i < 3; i++) {
    //    printf("%lf\t", vtr1[i]);
    //}
    //printf("\n");
    //for (int i = 0; i < 3; i++) {
    //    printf("%lf\t", vtr2[i]);
    //}
    //printf("\n");
    //

    long int dims[] = {1, 3};
    PyObject* resVtr = PyArray_SimpleNew(2, (npy_intp*)dims, NPY_DOUBLE);
    npy_intp* np_strides = PyArray_STRIDES(resVtr);
    char* pResVtr = (char*)PyArray_DATA(resVtr);
    
    *(double*)(pResVtr + 0 * np_strides[1]) = vtr1[1] * vtr2[2] - vtr1[2] * vtr2[1];
    *(double*)(pResVtr + 1 * np_strides[1]) = vtr1[2] * vtr2[0] - vtr1[0] * vtr2[2];
    *(double*)(pResVtr + 2 * np_strides[1]) = vtr1[0] * vtr2[1] - vtr1[1] * vtr2[0];
    
    //for (int i = 0; i < 3; i++) {
    //    printf("%lf\t", *(double*)(pResVtr + i * np_strides[1]));
    //}
    //printf("\n");

    return resVtr;
}

PyObject* VTR_ADD_3D(pVtr1, strd1, pVtr2, strd2) {
    double vtr1[] = {*(double*)((pVtr1) + (strd1) * 0), *(double*)((pVtr1) + (strd1) * 1), *(double*)((pVtr1) + (strd1) * 2)};
    double vtr2[] = {*(double*)((pVtr2) + (strd2) * 0), *(double*)((pVtr2) + (strd2) * 1), *(double*)((pVtr2) + (strd2) * 2)};

    long int dims[] = {1, 3};
    PyObject* resVtr = PyArray_SimpleNew(2, (npy_intp*)dims, NPY_DOUBLE);
    npy_intp* np_strides = PyArray_STRIDES(resVtr);
    char* pResVtr = (char*)PyArray_DATA(resVtr);

    *(double*)(pResVtr + 0 * np_strides[1]) = *(double*)((pVtr1) + (strd1) * 0) + *(double*)((pVtr2) + (strd2) * 0);
    *(double*)(pResVtr + 1 * np_strides[1]) = *(double*)((pVtr1) + (strd1) * 1) + *(double*)((pVtr2) + (strd2) * 1);
    *(double*)(pResVtr + 2 * np_strides[1]) = *(double*)((pVtr1) + (strd1) * 2) + *(double*)((pVtr2) + (strd2) * 2);
    
    return resVtr;
}

PyObject* VTR_SUBTRACT_3D(pVtr1, strd1, pVtr2, strd2) {
    double vtr1[] = {*(double*)((pVtr1) + (strd1) * 0), *(double*)((pVtr1) + (strd1) * 1), *(double*)((pVtr1) + (strd1) * 2)};
    double vtr2[] = {*(double*)((pVtr2) + (strd2) * 0), *(double*)((pVtr2) + (strd2) * 1), *(double*)((pVtr2) + (strd2) * 2)};

    long int dims[] = {1, 3};
    PyObject* resVtr = PyArray_SimpleNew(2, (npy_intp*)dims, NPY_DOUBLE);
    npy_intp* np_strides = PyArray_STRIDES(resVtr);
    char* pResVtr = (char*)PyArray_DATA(resVtr);

    *(double*)(pResVtr + 0 * np_strides[1]) = *(double*)((pVtr1) + (strd1) * 0) - *(double*)((pVtr2) + (strd2) * 0);
    *(double*)(pResVtr + 1 * np_strides[1]) = *(double*)((pVtr1) + (strd1) * 1) - *(double*)((pVtr2) + (strd2) * 1);
    *(double*)(pResVtr + 2 * np_strides[1]) = *(double*)((pVtr1) + (strd1) * 2) - *(double*)((pVtr2) + (strd2) * 2);
    
    return resVtr;
}


/*
    Real accelerations here!
*/
static PyObject* C_SurroundVectorSet(PyObject* self, PyObject* args) {
    double low_cutoff, high_cutoff;
    PyObject* animoAtoms = NULL;
    PyObject* animoAtomsArray = NULL;

    if (!PyArg_ParseTuple(args, "Odd", &animoAtoms, &low_cutoff, &high_cutoff))
        return NULL;

    animoAtomsArray = PyArray_FROM_OTF(animoAtoms, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    npy_intp* dims = PyArray_DIMS(animoAtomsArray);
    npy_intp* strides = PyArray_STRIDES(animoAtomsArray);
    char* pAnimoAtomsArray = (char*)PyArray_DATA(animoAtomsArray);

    // DEBUGGING
    printf("Cutoff = (%lf, %lf)\nDims:Strides = ", low_cutoff, high_cutoff); 
    int i;
    for (i = 0; i < PyArray_NDIM(animoAtomsArray); i++) {
        printf("%ld:%ld    ", dims[i], strides[i]);
    } 
    printf("\n");

    //i = 0;
    //int j, k;
    //for (j = 0; j < dims[1]; j++) {
    //    for (k = 0; k < dims[2]; k++) {
    //        //printf("%lf\t", *(double*)(pAnimoAtomsArray + i * strides[0] + j * strides[1] + k * strides[2]));
    //        printf("%lf\t", *PTR_DOUBLE_ELEM_3(pAnimoAtomsArray, strides, i, j, k));
    //    }
    //    printf("\n");
    //}

    //printf("%lf\n", VTR_DOT_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 0), strides[2], \
    //                           (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 1), strides[2]));

    //printf("%lf\n", VTR_LEN_3D((char*)pAnimoAtomsArray, strides[2]));
    //PyObject* crossRes = VTR_CROSS_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 0), strides[2], \
    //                                  (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 2), strides[2]);
    //for (i = 0; i < 3; i++) {
    //    printf("%lf\t", *PTR_DOUBLE_ELEM_2(PyArray_DATA(crossRes), PyArray_STRIDES(crossRes), 0, i));
    //}
    //printf("\n");

    //printf("Adding Vtrs: ");
    //PyObject* addRes = VTR_ADD_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 0), strides[2], \
    //                                  (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 2), strides[2]);
    //for (i = 0; i < 3; i++) {
    //    printf("%lf\t", *PTR_DOUBLE_ELEM_2(PyArray_DATA(addRes), PyArray_STRIDES(addRes), 0, i));
    //}
    //printf("\n");
    
    //printf("Subtracting Vtrs: ");
    //PyObject* subRes = VTR_SUBTRACT_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 0), strides[2], \
    //                                  (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, 0, 2), strides[2]);
    //for (i = 0; i < 3; i++) {
    //    printf("%lf\t", *PTR_DOUBLE_ELEM_2(PyArray_DATA(subRes), PyArray_STRIDES(subRes), 0, i));
    //}
    //printf("\n");

    // END OF DEBUGGING

    PyObject* aminosList = PyList_New(0);
    int aa1, aa2;
    for (aa1 = 0; aa1 < dims[0]; aa1++) {
        PyObject* theAminoContact = PyList_New(0);
        PyObject* vtr_CaN = VTR_SUBTRACT_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, aa1, 1), strides[2], \
                                            (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, aa1, 0), strides[2]);
        PyObject* vtr_CaC = VTR_SUBTRACT_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, aa1, 1), strides[2], \
                                            (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, aa1, 2), strides[2]);

        //if (aa1 == 0) {
        //    printf("CaN, CaC: \n");
        //    PRINT_VTR((char*)PyArray_DATA(vtr_CaN), PyArray_STRIDE(vtr_CaN, 1));
        //    PRINT_VTR((char*)PyArray_DATA(vtr_CaC), PyArray_STRIDE(vtr_CaC, 1));
        //    printf("\n");
        //}

        for (aa2 = 0; aa2 < dims[0]; aa2++) {
            PyObject* vtr_CaCa = VTR_SUBTRACT_3D((char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, aa2, 1), strides[2], \
                                                 (char*)PTR_DOUBLE_ELEM_2(pAnimoAtomsArray, strides, aa1, 1), strides[2]);
            double theVtrLen = VTR_LEN_3D((char*)PyArray_DATA(vtr_CaCa), PyArray_STRIDE(vtr_CaCa, 1));

            //if (aa1 == 0 && aa2 == 3) {
            //    printf("CaCa: \n");
            //    PRINT_VTR((char*)PyArray_DATA(vtr_CaCa), PyArray_STRIDE(vtr_CaCa, 1));
            //    printf("theVtrLen=%lf\n", theVtrLen);
            //}

            if (theVtrLen > low_cutoff && theVtrLen < high_cutoff) {
                double bondAngle = acos(VTR_ANGLECOS((char*)PyArray_DATA(vtr_CaN), PyArray_STRIDE(vtr_CaN, 1), \
                                                //(char*)PyArray_DATA(vtr_CaCa), PyArray_STRIDE(vtr_CaCa, 1), (aa1 == 0 && aa2 == 3));
                                                (char*)PyArray_DATA(vtr_CaCa), PyArray_STRIDE(vtr_CaCa, 1)));
                
                if (aa1 == 0 && aa2 == 3) {
                    printf("ANGLECOS = %lf\n", bondAngle);
                    printf("STRIDE = %ld\n", PyArray_STRIDE(vtr_CaN, 1));
                }

                PyObject* plnNorm1 = VTR_CROSS_3D( \
                    (char*)PyArray_DATA(vtr_CaN), PyArray_STRIDE(vtr_CaN, 1), \
                    (char*)PyArray_DATA(vtr_CaCa), PyArray_STRIDE(vtr_CaCa, 1)  \
                );
                PyObject* plnNorm2 = VTR_CROSS_3D( \
                    (char*)PyArray_DATA(vtr_CaN), PyArray_STRIDE(vtr_CaN, 1), \
                    (char*)PyArray_DATA(vtr_CaC), PyArray_STRIDE(vtr_CaC, 1)  \
                );
                double torsionAngle = acos(VTR_ANGLECOS((char*)PyArray_DATA(plnNorm1), PyArray_STRIDE(plnNorm1, 1), \
                                                   (char*)PyArray_DATA(plnNorm2), PyArray_STRIDE(plnNorm2, 1)));

                PyObject* contactInfo = PyList_New(3);
                PyList_SetItem(contactInfo, 0, PyFloat_FromDouble(theVtrLen));
                PyList_SetItem(contactInfo, 1, PyFloat_FromDouble(bondAngle));
                PyList_SetItem(contactInfo, 2, PyFloat_FromDouble(torsionAngle));
                PyList_Append(theAminoContact, contactInfo);
            }
        }
        PyList_Append(aminosList, theAminoContact);
    }


    //Py_RETURN_NONE;
    return aminosList;
}
    

static PyMethodDef ContactAccelMethods[] = {
    {"C_SurroundVectorSet", C_SurroundVectorSet, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef contactAccel_module = {
    PyModuleDef_HEAD_INIT,
    "ContactAccel", 
    NULL, 
    -1, 
    ContactAccelMethods
};

int init_numpy() {
    import_array();
    return 0;
}

PyMODINIT_FUNC
PyInit_ContactAccel(void) {
    init_numpy();
    return PyModule_Create(&contactAccel_module);
}
