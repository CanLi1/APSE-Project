#include "wrapper.h"
#include "stdio.h"
PyObject *make_int_list(int array[], size_t size) {

    PyObject *l = PyList_New(size);
    for (size_t i = 0; i != size; ++i) 
    {
        PyList_SET_ITEM(l, i, PyInt_FromLong(array[i]));
    }
    return l;

}

PyObject *make_real_list(real array[], size_t size) {

    PyObject *l = PyList_New(size);
    for (size_t i = 0; i != size; ++i) 
    {
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(array[i]));
    }
    return l;
}

void cleanup_array(PyObject *tup, int size) 
{
    for (int i = 1 ; i <= size ; i++) 
    {
        Py_DECREF(PyTuple_GetItem(tup, i));
    }
}

void pargs_dec_ref(PyObject *tup) {
    // Py_DECREF(PyTuple_GetItem(tup, 0));
    Py_DECREF(PyTuple_GetItem(tup, 1));
    Py_DECREF(PyTuple_GetItem(tup, 2));
    Py_DECREF(PyTuple_GetItem(tup, 3));
    Py_DECREF(PyTuple_GetItem(tup, 4));

}

real py_wrapper(arglist *al, char* module_name , char* funcname)
{
  
    /*
     *  Ideal gas molar volume, and AD example
     * 
     * */
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    PyObject *pValue_derivs , *pValue_hes;

    if (!Py_IsInitialized()) 
    {
        Py_Initialize();
    }

    pName = PyString_FromString(module_name);

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    //Initialize return_value.
    static real return_value = 0 ;

    int total_args = al->n; // length of args in at.
    int hessian_size = ((al->n * al->n) - al->n)/2 + al->n; //(n*n-n)/2 + n.
    
    if (al->at) total_args += 1;
    if (al->hes) total_args += 1;
    if (al->derivs) total_args += 1;
    
    if (pModule != NULL) 
    {
        pFunc = PyObject_GetAttrString(pModule,funcname);
        if (pFunc && PyCallable_Check(pFunc)) 
        {
            pArgs = PyTuple_New(total_args);
            // 3 arguments, always added: n, at, ra.
            PyTuple_SetItem(pArgs, 0, PyInt_FromLong(al->n));
            PyTuple_SetItem(pArgs, 1, make_int_list(al->at, al->n));
            PyTuple_SetItem(pArgs, 2, make_real_list(al->ra,al->n)); // Is al->ra always of length n?
            
            // 2 optional lists: derivatives , Hessian matrix.
            if (al->derivs) PyTuple_SetItem(pArgs, 3, make_real_list(al->derivs, al->n));
            if (al->hes) PyTuple_SetItem(pArgs, 4, make_real_list(al->hes,hessian_size));

            pValue = PyObject_CallObject(pFunc, pArgs);

            if (pValue != NULL) 
            {
                if (PyList_Check(pValue)) 
                {
                    // Returned values : V, derivs , hes
                    // Item 0 : V
                    return_value = PyFloat_AsDouble(PyList_GetItem(pValue,0));
       
                    // Item 1 (optional) : derivs
                    if (al->derivs) 
                    {
                        if (PyList_Check(PyList_GetItem(pValue,1))) 
                        {
                            //pValue_derivs contains the list we are looking for:
                            pValue_derivs = PyList_GetItem(pValue,1);
                            for (int i = 0 ; i < al->n ; i++) 
                            {
                                al->derivs[i] = PyFloat_AsDouble(PyList_GetItem(pValue_derivs,i));
                            }
                            Py_DECREF(pValue_derivs);
                        }
                    }
                    
                    // Item 2 (optional) : Hessian matrix
                    if (al->hes) 
                    {
                        if (PyList_Check(PyList_GetItem(pValue,2))) 
                        {
                            //pValue_hes contains the list we are looking for:
                            pValue_hes = PyList_GetItem(pValue,2);
                            for (int i = 0 ; i < hessian_size ; i++) 
                            {
                                al->hes[i] = PyFloat_AsDouble(PyList_GetItem(pValue_hes,i));
                            }
                            Py_DECREF(pValue_hes);
                        }
                    }
                }
                Py_DECREF(pValue);
            }

            else 
            {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                pargs_dec_ref(pArgs);
                Py_XDECREF(pArgs);
                PyErr_Print();
                return BAD_RETURN_VALUE;
            }
        }
        else 
        {
            return BAD_FUNCTION_DEFINITION;
        }
    }
    else 
    {
        return BAD_MODULE_NAME;
    }
    return return_value;

}