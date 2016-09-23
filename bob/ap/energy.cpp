/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Thu  6 Feb 09:00:05 2014
 *
 * @brief Bindings to the base class bob::ap::Energy
 */

#include <bob.blitz/cppapi.h>
#include <bob.blitz/cleanup.h>
#include <bob.extension/defines.h>
#include "types.h"

PyDoc_STRVAR(s_energy_str, BOB_EXT_MODULE_PREFIX ".Energy");

PyDoc_STRVAR(s_energy_doc,
"Energy(sampling_frequency, [win_length_ms=20., [win_shift_ms=10., [normalize_mean=True]]]) -> new Energy\n\
Energy(other) -> new Energy\n\
\n\
Objects of this class, after configuration, can extract the energy\n\
of frames extracted from a 1D audio array/signal.\n\
\n\
Parameters:\n\
\n\
sampling_frequency\n\
  [float] the sampling frequency/frequency rate\n\
\n\
win_length_ms\n\
  [float] the window length in miliseconds\n\
\n\
win_shift_ms\n\
  [float] the window shift in miliseconds\n\
\n\
normalize_mean\n\
  [bool] Tells whether frame should be normalized \n\
  by subtracting mean (True) or dividing by max_range (False)\n\
  ``True`` is the default value.\n\
\n\
other\n\
  [Energy] an object of which is or inherits from ``Energy``\n\
  that will be deep-copied into a new instance.\n\
\n\
"
);

int PyBobApEnergy_Check(PyObject* o) {
  return PyObject_IsInstance(o, reinterpret_cast<PyObject*>(&PyBobApEnergy_Type));
}

static void PyBobApEnergy_Delete (PyBobApEnergyObject* o) {

  o->parent.cxx = 0;
  delete o->cxx;
  Py_TYPE(o)->tp_free((PyObject*)o);

}

static int PyBobApEnergy_InitCopy
(PyBobApEnergyObject* self, PyObject* args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"other", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  PyObject* other = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
        &PyBobApEnergy_Type, &other)) return -1;

  auto copy = reinterpret_cast<PyBobApEnergyObject*>(other);

  try {
    self->cxx = new bob::ap::Energy(*(copy->cxx));
    if (!self->cxx) {
      PyErr_Format(PyExc_MemoryError, "cannot create new object of type `%s' - no more memory", Py_TYPE(self)->tp_name);
      return -1;
    }
    self->parent.cxx = self->cxx;
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot create new object of type `%s' - unknown exception thrown", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

static int PyBobApEnergy_InitParameters
(PyBobApEnergyObject* self, PyObject *args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {
    "sampling_frequency",
    "win_length_ms",
    "win_shift_ms",
    "normalize_mean",
    0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  double sampling_frequency = 0.;
  double win_length_ms = 20.;
  double win_shift_ms = 10.;
  PyObject* normalize_mean = Py_True;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|ddO", kwlist,
        &sampling_frequency, &win_length_ms, &win_shift_ms, &normalize_mean)) return -1;

  bool normalize_mean_ = PyObject_IsTrue(normalize_mean);

  try {
    self->cxx = new bob::ap::Energy(sampling_frequency,
        win_length_ms, win_shift_ms, normalize_mean_);
    if (!self->cxx) {
      PyErr_Format(PyExc_MemoryError, "cannot create new object of type `%s' - no more memory", Py_TYPE(self)->tp_name);
      return -1;
    }
    self->parent.cxx = self->cxx;
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot create new object of type `%s' - unknown exception thrown", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0; ///< SUCCESS

}

static int PyBobApEnergy_Init(PyBobApEnergyObject* self,
    PyObject* args, PyObject* kwds) {

  Py_ssize_t nargs = (args?PyTuple_Size(args):0) + (kwds?PyDict_Size(kwds):0);

  switch (nargs) {

    case 1:

      {

        PyObject* arg = 0; ///< borrowed (don't delete)
        if (PyTuple_Size(args)) arg = PyTuple_GET_ITEM(args, 0);
        else {
          PyObject* tmp = PyDict_Values(kwds);
          auto tmp_ = make_safe(tmp);
          arg = PyList_GET_ITEM(tmp, 0);
        }

        if (PyBobApEnergy_Check(arg)) {
          return PyBobApEnergy_InitCopy(self, args, kwds);
        }

        else {
          return PyBobApEnergy_InitParameters(self, args, kwds);
        }

        PyErr_Format(PyExc_TypeError, "cannot initialize `%s' with `%s' (see help)", Py_TYPE(self)->tp_name, Py_TYPE(arg)->tp_name);

      }

      break;

    default:

      return PyBobApEnergy_InitParameters(self, args, kwds);

  }

  return -1;

}

static PyObject* PyBobApEnergy_Repr(PyBobApEnergyObject* self) {
  static const int MAXSIZE = 256;
  char buffer[MAXSIZE];
  auto count = std::snprintf(buffer, MAXSIZE, "%s(sampling_frequency=%f, win_length_ms=%f, win_shift_ms=%f, normalize_mean=%s)", Py_TYPE(self)->tp_name, self->cxx->getSamplingFrequency(), self->cxx->getWinLengthMs(), self->cxx->getWinShiftMs(), self->cxx->getNormalizeMean()?"True":"False");
  return
# if PY_VERSION_HEX >= 0x03000000
  PyUnicode_FromStringAndSize
# else
  PyString_FromStringAndSize
# endif
    (buffer, (count<=MAXSIZE)?count:MAXSIZE);
}

static PyObject* PyBobApEnergy_RichCompare (PyBobApEnergyObject* self,
    PyObject* other, int op) {

  if (!PyBobApEnergy_Check(other)) {
    PyErr_Format(PyExc_TypeError, "cannot compare `%s' with `%s'",
        Py_TYPE(self)->tp_name, Py_TYPE(other)->tp_name);
    return 0;
  }

  auto other_ = reinterpret_cast<PyBobApEnergyObject*>(other);

  switch (op) {
    case Py_EQ:
      if (self->cxx->operator==(*other_->cxx)) Py_RETURN_TRUE;
      Py_RETURN_FALSE;
      break;
    case Py_NE:
      if (self->cxx->operator!=(*other_->cxx)) Py_RETURN_TRUE;
      Py_RETURN_FALSE;
      break;
    default:
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
  }

}

PyDoc_STRVAR(s_energy_floor_str, "energy_floor");
PyDoc_STRVAR(s_energy_floor_doc,
"The energy flooring threshold"
);

static PyObject* PyBobApEnergy_GetEnergyFloor
(PyBobApEnergyObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getEnergyFloor());
}

static int PyBobApEnergy_SetEnergyFloor
(PyBobApEnergyObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' energy floor can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setEnergyFloor(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `energy_floor' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

static PyGetSetDef PyBobApEnergy_getseters[] = {
    {
      s_energy_floor_str,
      (getter)PyBobApEnergy_GetEnergyFloor,
      (setter)PyBobApEnergy_SetEnergyFloor,
      s_energy_floor_doc,
      0
    },
    {0}  /* Sentinel */
};

static PyObject* PyBobApEnergy_Call
(PyBobApEnergyObject* self, PyObject *args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"input", "output", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  PyBlitzArrayObject* input = 0;
  PyBlitzArrayObject* output = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O&|O&", kwlist,
        &PyBlitzArray_Converter, &input,
        &PyBlitzArray_OutputConverter, &output
        )) return 0;

  auto input_ = make_safe(input);
  auto output_ = make_xsafe(output);

  if (input->type_num != NPY_FLOAT64) {
    PyErr_Format(PyExc_TypeError, "`%s' only supports 1D 64-bit float arrays for input array `input'", Py_TYPE(self)->tp_name);
    return 0;
  }

  if (input->ndim != 1) {
    PyErr_Format(PyExc_TypeError, "`%s' only supports 1D 64-bit float arrays for input array `input'", Py_TYPE(self)->tp_name);
    return 0;
  }

  auto bz_input = PyBlitzArrayCxx_AsBlitz<double,1>(input);

  if (output) {

    if (output->type_num != NPY_FLOAT64) {
      PyErr_Format(PyExc_TypeError, "`%s' only supports 64-bit float arrays for output array `output'", Py_TYPE(self)->tp_name);
      return 0;
    }

    if (input->ndim != output->ndim) {
      PyErr_Format(PyExc_RuntimeError, "Input and output arrays should have matching number of dimensions, but input array `input' has %" PY_FORMAT_SIZE_T "d dimensions while output array `output' has %" PY_FORMAT_SIZE_T "d dimensions", input->ndim, output->ndim);
      return 0;
    }

  }

  else {

    Py_ssize_t length = self->cxx->getShape(*bz_input)(0);
    output = (PyBlitzArrayObject*)PyBlitzArray_SimpleNew(NPY_FLOAT64, 1, &length);
    if (!output) return 0;
    output_ = make_safe(output);

  }

  auto bz_output = PyBlitzArrayCxx_AsBlitz<double,1>(output);

  try {
    self->cxx->operator()(*bz_input, *bz_output);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return 0;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot call object of type `%s' - unknown exception thrown", Py_TYPE(self)->tp_name);
    return 0;
  }

  return PyBlitzArray_NUMPY_WRAP(Py_BuildValue("O", output));

}

PyTypeObject PyBobApEnergy_Type = {
    PyVarObject_HEAD_INIT(0, 0)
    s_energy_str,                             /*tp_name*/
    sizeof(PyBobApEnergyObject),              /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)PyBobApEnergy_Delete,         /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    (reprfunc)PyBobApEnergy_Repr,             /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash */
    (ternaryfunc)PyBobApEnergy_Call,          /* tp_call */
    (reprfunc)PyBobApEnergy_Repr,             /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    s_energy_doc,                             /* tp_doc */
    0,		                                    /* tp_traverse */
    0,		                                    /* tp_clear */
    (richcmpfunc)PyBobApEnergy_RichCompare,    /* tp_richcompare */
    0,		                                    /* tp_weaklistoffset */
    0,		                                    /* tp_iter */
    0,		                                    /* tp_iternext */
    0,                                        /* tp_methods */
    0,                                        /* tp_members */
    PyBobApEnergy_getseters,                  /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    (initproc)PyBobApEnergy_Init,             /* tp_init */
    0,                                        /* tp_alloc */
    0,                                        /* tp_new */
};
