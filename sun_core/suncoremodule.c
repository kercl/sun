/* Copyright (c) 2018 Clemens Kerschbaumer
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <stdlib.h>
#include <limits.h>
#include <python3.7m/Python.h>

#include "sun_core/irrep.h"

int
read_dynkin_int(PyObject *obj, gt_int_t *v) {
    if (!PyLong_Check(obj))
        return 0;

    long x = PyLong_AsLong(obj);
    if (x > GT_INT_MAX || x < 0)
        return 0;

    *v = (gt_int_t)x;

    return 1;
}

static PyObject *
sun_core_irrep_dimension_from_dynkin(PyObject *self, PyObject *args) {
    PyObject *obj, *seq, *item;
    int len = 0;
    gt_int_t *dynkin = NULL;

    if (!PyArg_ParseTuple(args, "O", &obj))
        return Py_None;

    seq = PySequence_Fast(obj, "expected a sequence");
    len = PySequence_Size(obj);

    char value_error_msg[200];
    snprintf(value_error_msg, sizeof(value_error_msg),
        "Dynkin labels are supposed to be lists"\
        "of positive integers limited by [%d, %d].",
        GT_INT_MIN, GT_INT_MAX);

    if (PyList_Check(seq)) {
        dynkin = malloc(sizeof(gt_int_t) * len);

        for (int i = 0; i < len; i++) {
            item = PyList_GET_ITEM(seq, i);
            if (!read_dynkin_int(item, dynkin + i)) {
                free(dynkin);
                PyErr_SetString(PyExc_ValueError, value_error_msg);
                return NULL;
            }
        }
    } else {
        PyErr_SetString(PyExc_ValueError, value_error_msg);
        return NULL;
    }
    Py_DECREF(seq);

    return Py_BuildValue("i", dimension_from_dynkin(dynkin, len));
}

static PyMethodDef sun_core_methods[] = {
    {"irrep_dimension_from_dynkin",
     sun_core_irrep_dimension_from_dynkin,
     METH_VARARGS,
     "For a given dynkin label, return the"\
     "dimension of the corresponding irreducible"\
     "representation."},

    {NULL,
     NULL,
     0,
     NULL}
};

static struct PyModuleDef sun_core = {
    PyModuleDef_HEAD_INIT,
    "sun_core",
    "This module provides the core functionality"\
    "for sun.",
    -1,
    sun_core_methods
};

PyMODINIT_FUNC PyInit_sun_core(void) {
    return PyModule_Create(&sun_core);
}
