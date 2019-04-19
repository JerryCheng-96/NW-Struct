#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* sys_cmd(PyObject* self, PyObject* args) {
    const char* command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return PyLong_FromLong(sts);
}

int test_func(int number) {
    return number + 2;
}

static PyObject* _test_func(PyObject* self, PyObject* args) {
    int number;

    if (!PyArg_ParseTuple(args, "i", &number))
        return NULL;

    int result = test_func(number);

    return PyLong_FromLong(result);
}

static PyObject* test_array(PyObject* self, PyObject* args) {
    PyObject *arg1=NULL, *arg2=NULL, *out=NULL;
    PyObject *arr1=NULL, *arr2=NULL, *oarr=NULL;

    if (!PyArg_ParseTuple(args, "OOO!", &arg1, &arg2,
        &PyArray_Type, &out)) return NULL;

    arr1 = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr1 == NULL) return NULL;
    arr2 = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
    if (arr2 == NULL) goto fail;
    oarr = PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (oarr == NULL) goto fail;

    /* code that makes use of arguments */
    /* You will probably need at least
       nd = PyArray_NDIM(<..>)    -- number of dimensions
       dims = PyArray_DIMS(<..>)  -- npy_intp array of length nd
                                     showing length in each dim.
       dptr = (double *)PyArray_DATA(<..>) -- pointer to data.

       If an error occurs goto fail.
     */

    //int nd = PyArray_NDIM(arr1);
    //npy_intp* np_dim = PyArray_DIMS(arr1);
    //npy_intp* np_strides = PyArray_STRIDES(arr1);
    //long int dims[] = {np_dim[0], np_dim[1]};
    //long int strides[] = {np_strides[0], np_strides[1]};

    //char* arr1_pData = (char*)PyArray_DATA(arr1);
    //for (int i; i < dims[0]; i++) {
    //    for (int j; j < dims[1]; j++) {
    //        printf("(%d, %d, %d, %d)=", i, j, dims[0], dims[1]);
    //        printf("%lf\t", *(double*)(arr1_pData + i * strides[0] + j * strides[1]));
    //    }
    //    printf("\n");
    //}
    
    //int nd = PyArray_NDIM(oarr);
    //npy_intp* np_dim = PyArray_DIMS(oarr);
    //npy_intp* np_strides = PyArray_STRIDES(oarr);
    //printf("%ld, %ld\n", np_dim[0], np_dim[1]);
    //char* oarr_pData = (char*)PyArray_DATA(oarr);
    //*(double*)(oarr_pData + np_strides[0]) = 0;

    long int dims[] = {2, 3};
    PyObject* newArr = PyArray_SimpleNew(2, (npy_intp*)dims, NPY_DOUBLE);
    npy_intp* np_strides = PyArray_STRIDES(newArr);
    char* pNewArr = (char*)PyArray_DATA(newArr);
    int i,j;
    double val = 1.0;
    for (i = 0; i < 2; i++){
        for (j = 0; j < 3; j++) {
            *(double*)(pNewArr + i * np_strides[0] + j * np_strides[1]) = val;
            val += 1;
        }
    }

    return newArr;

 fail:
    return NULL;
}

static PyMethodDef TestModMethods[] = {
    {"sys_cmd", sys_cmd, METH_VARARGS, ""},
    {"test_func", _test_func, METH_VARARGS, ""},
    {"test_array", test_array, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef test_module = {
    PyModuleDef_HEAD_INIT,
    "TestMod", 
    NULL, 
    -1, 
    TestModMethods
};

int init_numpy() {
    import_array();
    return 0;
}

PyMODINIT_FUNC
PyInit_TestMod(void) {
    init_numpy();
    return PyModule_Create(&test_module);
}
