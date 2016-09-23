/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Thu  6 Feb 09:00:05 2014
 *
 * @brief Bindings to the base class bob::ap::FrameExtractor
 */

#include <bob.blitz/cppapi.h>
#include <bob.blitz/cleanup.h>
#include <bob.extension/defines.h>
#include "types.h"

PyDoc_STRVAR(s_frame_extractor_str, BOB_EXT_MODULE_PREFIX ".FrameExtractor");

PyDoc_STRVAR(s_frame_extractor_doc,
"FrameExtractor(sampling_frequency, [win_length_ms=20., [win_shift_ms=10., [normalize_mean=True]]]) -> new FrameExtractor\n\
FrameExtractor(other) -> new FrameExtractor\n\
\n\
This class is a base type for classes that perform audio\n\
processing on a frame basis. It *can* be instantiated from Python.\n\
\n\
Objects of this class, after configuration, can extract audio frame\n\
from a 1D audio array/signal. You can instantiate objects of this\n\
class by passing a set of construction parameters or another object\n\
of which the base type is ``FrameExtractor``.\n\
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
  [FrameExtractor] an object of which is or inherits from a FrameExtractor\n\
  that will be deep-copied into a new instance.\n\
\n\
"
);

int PyBobApFrameExtractor_Check(PyObject* o) {
  return PyObject_IsInstance(o, reinterpret_cast<PyObject*>(&PyBobApFrameExtractor_Type));
}

static void PyBobApFrameExtractor_Delete (PyBobApFrameExtractorObject* o) {

  delete o->cxx;
  Py_TYPE(o)->tp_free((PyObject*)o);

}

static int PyBobApFrameExtractor_InitCopy
(PyBobApFrameExtractorObject* self, PyObject* args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"other", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  PyObject* other = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
        &PyBobApFrameExtractor_Type, &other)) return -1;

  auto copy = reinterpret_cast<PyBobApFrameExtractorObject*>(other);

  try {
    self->cxx = new bob::ap::FrameExtractor(*(copy->cxx));
    if (!self->cxx) {
      PyErr_Format(PyExc_MemoryError, "cannot create new object of type `%s' - no more memory", Py_TYPE(self)->tp_name);
      return -1;
    }
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

static int PyBobApFrameExtractor_InitParameters
(PyBobApFrameExtractorObject* self, PyObject *args, PyObject* kwds) {

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
    self->cxx = new bob::ap::FrameExtractor(sampling_frequency, 
        win_length_ms, win_shift_ms, normalize_mean_);
    if (!self->cxx) {
      PyErr_Format(PyExc_MemoryError, "cannot create new object of type `%s' - no more memory", Py_TYPE(self)->tp_name);
      return -1;
    }
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

static int PyBobApFrameExtractor_Init(PyBobApFrameExtractorObject* self,
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

        if (PyBobApFrameExtractor_Check(arg)) {
          return PyBobApFrameExtractor_InitCopy(self, args, kwds);
        }

        else {
          return PyBobApFrameExtractor_InitParameters(self, args, kwds);
        }

        PyErr_Format(PyExc_TypeError, "cannot initialize `%s' with `%s' (see help)", Py_TYPE(self)->tp_name, Py_TYPE(arg)->tp_name);

      }

      break;

    default:

      return PyBobApFrameExtractor_InitParameters(self, args, kwds);

  }

  return -1;

}

static PyObject* PyBobApFrameExtractor_Repr(PyBobApFrameExtractorObject* self) {
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

static PyObject* PyBobApFrameExtractor_RichCompare (PyBobApFrameExtractorObject* self,
    PyObject* other, int op) {

  if (!PyBobApFrameExtractor_Check(other)) {
    PyErr_Format(PyExc_TypeError, "cannot compare `%s' with `%s'",
        Py_TYPE(self)->tp_name, Py_TYPE(other)->tp_name);
    return 0;
  }

  auto other_ = reinterpret_cast<PyBobApFrameExtractorObject*>(other);

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

PyDoc_STRVAR(s_sampling_frequency_str, "sampling_frequency");
PyDoc_STRVAR(s_sampling_frequency_doc,
"The sampling frequency/frequency rate"
);

static PyObject* PyBobApFrameExtractor_GetSamplingFrequency
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getSamplingFrequency());
}

static int PyBobApFrameExtractor_SetSamplingFrequency
(PyBobApFrameExtractorObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' sampling frequency can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setSamplingFrequency(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `sampling_frequency' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_win_length_ms_str, "win_length_ms");
PyDoc_STRVAR(s_win_length_ms_doc,
"The window length of the cepstral analysis in milliseconds"
);

static PyObject* PyBobApFrameExtractor_GetWinLengthMs
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getWinLengthMs());
}

static int PyBobApFrameExtractor_SetWinLengthMs
(PyBobApFrameExtractorObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' windows length can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setWinLengthMs(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `win_length_ms' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_win_shift_ms_str, "win_shift_ms");
PyDoc_STRVAR(s_win_shift_ms_doc,
"The window shift of the cepstral analysis in milliseconds"
);

static PyObject* PyBobApFrameExtractor_GetWinShiftMs
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getWinShiftMs());
}

static int PyBobApFrameExtractor_SetWinShiftMs
(PyBobApFrameExtractorObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' windows shift can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setWinShiftMs(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `win_shift_ms' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_normalize_mean_str, "normalize_mean");
PyDoc_STRVAR(s_normalize_mean_doc,
"Tells whether frame should be normalized by subtracting mean (True) or dividing by max_range (False)\n\
");

static PyObject* PyBobApFrameExtractor_GetNormalizeMean
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  if (self->cxx->getNormalizeMean()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApFrameExtractor_SetNormalizeMean
(PyBobApFrameExtractorObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setNormalizeMean(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `normalize_mean' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;
}

PyDoc_STRVAR(s_win_length_str, "win_length");
PyDoc_STRVAR(s_win_length_doc,
"The normalized window length w.r.t. the sample frequency"
);

static PyObject* PyBobApFrameExtractor_GetWinLength
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  return Py_BuildValue("n", self->cxx->getWinLength());
}

PyDoc_STRVAR(s_win_shift_str, "win_shift");
PyDoc_STRVAR(s_win_shift_doc,
"The normalized window shift w.r.t. the sample frequency"
);

static PyObject* PyBobApFrameExtractor_GetWinShift
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  return Py_BuildValue("n", self->cxx->getWinShift());
}

static PyGetSetDef PyBobApFrameExtractor_getseters[] = {
    {
      s_sampling_frequency_str,
      (getter)PyBobApFrameExtractor_GetSamplingFrequency,
      (setter)PyBobApFrameExtractor_SetSamplingFrequency,
      s_sampling_frequency_doc,
      0
    },
    {
      s_win_length_ms_str,
      (getter)PyBobApFrameExtractor_GetWinLengthMs,
      (setter)PyBobApFrameExtractor_SetWinLengthMs,
      s_win_length_ms_doc,
      0
    },
    {
      s_win_shift_ms_str,
      (getter)PyBobApFrameExtractor_GetWinShiftMs,
      (setter)PyBobApFrameExtractor_SetWinShiftMs,
      s_win_shift_ms_doc,
      0
    },
    {
      s_win_length_str,
      (getter)PyBobApFrameExtractor_GetWinLength,
      0,
      s_win_length_doc,
      0
    },
    {
      s_win_shift_str,
      (getter)PyBobApFrameExtractor_GetWinShift,
      0,
      s_win_shift_doc,
      0
    },
    {
      s_normalize_mean_str,
      (getter)PyBobApFrameExtractor_GetNormalizeMean,
      (setter)PyBobApFrameExtractor_SetNormalizeMean,
      s_normalize_mean_doc,
      0
    },
    {0}  /* Sentinel */
};

PyDoc_STRVAR(s_shape_str, "get_shape");
PyDoc_STRVAR(s_shape_doc,
"x.get_shape(input) -> tuple\n\
\n\
Computes the shape of the output features, given the size of\n\
an input array or an input array.\n\
\n\
Parameters:\n\
\n\
input\n\
  [int|array] Either an integral value or an array for which\n\
  the output shape of this extractor is going to be computed.\n\
\n\
This method always returns a 2-tuple containing the shape of\n\
output features produced by this extractor.\n\
");

static PyObject* PyBobApFrameExtractor_GetShapeInt
(PyBobApFrameExtractorObject* self, PyObject* args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"input", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  Py_ssize_t input = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "n", kwlist, &input)) return 0;

  blitz::TinyVector<int,2> retval = self->cxx->getShape(input);

  return Py_BuildValue("(nn)", retval[0], retval[1]);

}

static PyObject* PyBobApFrameExtractor_GetShapeArray
(PyBobApFrameExtractorObject* self, PyObject* args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"input", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  PyBlitzArrayObject* input = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O&", kwlist,
        &input, &PyBlitzArray_Converter)) return 0;
  auto input_ = make_safe(input);

  if (input->ndim != 1 || input->type_num != NPY_FLOAT64) {
    PyErr_Format(PyExc_TypeError, "`%s' only accepts 1-dimensional 64-bit float arrays (not %" PY_FORMAT_SIZE_T "dD %s arrays)", Py_TYPE(self)->tp_name, input->ndim, PyBlitzArray_TypenumAsString(input->type_num));
    return 0;
  }

  blitz::TinyVector<int,2> retval =
    self->cxx->getShape(*PyBlitzArrayCxx_AsBlitz<double,1>(input));

  return Py_BuildValue("(nn)", retval[0], retval[1]);

}

static PyObject* PyBobApFrameExtractor_GetShape
(PyBobApFrameExtractorObject* self, PyObject* args, PyObject* kwds) {

  // input object can either be an integral value or an array<double,1>
  // the return value is always a 2-tuple

  Py_ssize_t nargs = (args?PyTuple_Size(args):0) + (kwds?PyDict_Size(kwds):0);

  if (nargs != 1) {
    PyErr_Format(PyExc_RuntimeError, "%s.%s expects 1 parameter, but you passed %" PY_FORMAT_SIZE_T "d", Py_TYPE(self)->tp_name, s_shape_str, nargs);
    return 0;
  }

  PyObject* arg = 0; ///< borrowed (don't delete)
  if (PyTuple_Size(args)) arg = PyTuple_GET_ITEM(args, 0);
  else {
    PyObject* tmp = PyDict_Values(kwds);
    auto tmp_ = make_safe(tmp);
    arg = PyList_GET_ITEM(tmp, 0);
  }

  if (PyInt_Check(arg)) {
    return PyBobApFrameExtractor_GetShapeInt(self, args, kwds);
  }

  return PyBobApFrameExtractor_GetShapeArray(self, args, kwds);

}

static PyMethodDef PyBobApFrameExtractor_methods[] = {
    {
      s_shape_str,
      (PyCFunction)PyBobApFrameExtractor_GetShape,
      METH_VARARGS|METH_KEYWORDS,
      s_shape_doc
    },
    {0}  /* Sentinel */
};

PyTypeObject PyBobApFrameExtractor_Type = {
    PyVarObject_HEAD_INIT(0, 0)
    s_frame_extractor_str,                    /*tp_name*/
    sizeof(PyBobApFrameExtractorObject),      /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)PyBobApFrameExtractor_Delete, /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    (reprfunc)PyBobApFrameExtractor_Repr,     /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash */
    0,                                        /* tp_call */
    (reprfunc)PyBobApFrameExtractor_Repr,     /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    s_frame_extractor_doc,                    /* tp_doc */
    0,		                                    /* tp_traverse */
    0,		                                    /* tp_clear */
    (richcmpfunc)PyBobApFrameExtractor_RichCompare,    /* tp_richcompare */
    0,		                                    /* tp_weaklistoffset */
    0,		                                    /* tp_iter */
    0,		                                    /* tp_iternext */
    PyBobApFrameExtractor_methods,            /* tp_methods */
    0,                                        /* tp_members */
    PyBobApFrameExtractor_getseters,          /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    (initproc)PyBobApFrameExtractor_Init,     /* tp_init */
    0,                                        /* tp_alloc */
    0,                                        /* tp_new */
};
