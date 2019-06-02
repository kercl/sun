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

#include <python3.7m/Python.h>

#include "sun_core/irrep.h"

static PyObject *
sun_core_irrep_dimension_from_dynkin(PyObject *self, PyObject *args) {
    printf("Hello World\n");
    return Py_None;
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

static struct PyModuleDef myModule = {
    PyModuleDef_HEAD_INIT,
    "sun_core",
    "This module provides the core functionality"\
    "for sun.",
    -1,
    sun_core_methods
};
