#include "Python.h"
#include <stdio.h>

#define No_AE_redefs
#include "solvers/funcadd.h"

// #define Printf printf

// Constants defining return values:
#define BAD_MODULE_NAME -1
#define BAD_RETURN_VALUE -2
#define BAD_FUNCTION_DEFINITION -3

// FIXME: Will not work without No_AE_redefs defined above since funcadd.h redefines getenv.
// #define PYTHON_MODULE_NAME getenv("WRAP_MOD")

#define str(s) #s

#define FUNCGEN(module_name,funcname) \
extern "C" real funcname (arglist *al) \
{  \
    return py_wrapper(al, str(module_name),str(funcname)); \
}

#define FUNCGEN_DUMMY(funcname) \
real funcname (arglist *al) \
{  \
    return py_wrapper_dummy(al, str(funcname)); \
}

#ifdef __cplusplus
extern "C" {
#endif
real py_wrapper(arglist *al, char* module_name , char* funcname);
real py_wrapper_dummy(arglist *al, char* funcname);
PyObject *make_int_list(int array[], size_t size);
PyObject *make_real_list(real array[], size_t size);
void cleanup_array(PyObject *tup, int size);
void pargs_dec_ref(PyObject *tup);
#ifdef __cplusplus
}
#endif