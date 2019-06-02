#include <python3.7m/Python.h>

static PyObject *ErrorObject;

typedef struct {
    PyObject_HEAD
    PyObject            *attr;        /* Attributes dictionary */
} SunObject;

static PyTypeObject Sun_Type;

#define SunObject_Check(v)      (Py_TYPE(v) == &Sun_Type)

static SunObject *
newSunObject(PyObject *arg)
{
    SunObject *self;
    self = PyObject_New(SunObject, &Sun_Type);
    if (self == NULL)
        return NULL;
    self->attr = NULL;
    return self;
}

/* Sun methods */

static void
Sun_dealloc(SunObject *self)
{
    Py_XDECREF(self->attr);
    PyObject_Del(self);
}

static PyObject *
Sun_demo(SunObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, ":demo"))
        return NULL;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef Sun_methods[] = {
    {"demo", (PyCFunction)Sun_demo,  METH_VARARGS,
        PyDoc_STR("demo() -> None")},
    {NULL, NULL}           /* sentinel */
};

static PyObject *
Sun_getattro(SunObject *self, PyObject *name)
{
    if (self->attr != NULL) {
        PyObject *v = PyDict_GetItemWithError(self->attr, name);
        if (v != NULL) {
            Py_INCREF(v);
            return v;
        }
        else if (PyErr_Occurred()) {
            return NULL;
        }
    }
    return PyObject_GenericGetAttr((PyObject *)self, name);
}

static int
Sun_setattr(SunObject *self, const char *name, PyObject *v)
{
    if (self->attr == NULL) {
        self->attr = PyDict_New();
        if (self->attr == NULL)
            return -1;
    }
    if (v == NULL) {
        int rv = PyDict_DelItemString(self->attr, name);
        if (rv < 0 && PyErr_ExceptionMatches(PyExc_KeyError))
            PyErr_SetString(PyExc_AttributeError,
                "delete non-existing Sun attribute");
        return rv;
    }
    else
        return PyDict_SetItemString(self->attr, name, v);
}

static PyTypeObject Sun_Type = {
    /* The ob_type field must be initialized in the module init function
     * to be portable to Windows without using C++. */
    PyVarObject_HEAD_INIT(NULL, 0)
    "sunmodule.Sun",            /*tp_name*/
    sizeof(SunObject),          /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)Sun_dealloc,    /*tp_dealloc*/
    0,                          /*tp_vectorcall_offset*/
    (getattrfunc)0,             /*tp_getattr*/
    (setattrfunc)Sun_setattr,   /*tp_setattr*/
    0,                          /*tp_as_async*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    (getattrofunc)Sun_getattro, /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    Sun_methods,                /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
};

static PyObject *
sun_from_gt_top_row(PyObject *self, PyObject *args)
{
    SunObject *rv;

    if (!PyArg_ParseTuple(args, ":new"))
        return NULL;
    rv = newSunObject(args);
    if (rv == NULL)
        return NULL;
    return (PyObject *)rv;
}

static PyMethodDef sun_methods[] = {
    {"from_gt_top_row", sun_from_gt_top_row, METH_VARARGS,
        PyDoc_STR("from_gt_top_row() -> generate su(n) representation from top row of Gelfant-Tsetlin pattern.")},
    {NULL,              NULL}           /* sentinel */
};

PyDoc_STRVAR(module_doc,
    "This is a template module just for instruction.");

static int
sun_exec(PyObject *m)
{
    /* Finalize the type object including setting type of the new type
     * object; doing it here is required for portability, too. */
    if (PyType_Ready(&Sun_Type) < 0)
        goto fail;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("xx.error", NULL, NULL);
        if (ErrorObject == NULL)
            goto fail;
    }
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);

    return 0;
 fail:
    Py_XDECREF(m);
    return -1;
}

static struct PyModuleDef_Slot sun_slots[] = {
    {Py_mod_exec, sun_exec},
    {0, NULL},
};

static struct PyModuleDef xxmodule = {
    PyModuleDef_HEAD_INIT,
    "xx",
    module_doc,
    0,
    sun_methods,
    sun_slots,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_sun(void)
{
    return PyModuleDef_Init(&xxmodule);
}
