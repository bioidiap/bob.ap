/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Thu  6 Feb 09:00:05 2014
 *
 * @brief Bindings to the base class bob::ap::FrameExtractor
 */

#include <xbob.blitz/cppapi.h>
#include <xbob.blitz/cleanup.h>
#include <bob/ap/FrameExtractor.h>

PyDoc_STRVAR(s_frame_extractor_str, XBOB_EXT_MODULE_PREFIX ".FrameExtractor");

PyDoc_STRVAR(s_frame_extractor_doc,
"FrameExtractor(sampling_frequency, [win_length_ms=20., [win_shift_ms=10.]]) -> new FrameExtractor\n\
FrameExtractor(other) -> new FrameExtractor\n\
\n\
This class is a base type for classes that perform audio processing on a\n\
frame basis. It *can* be instantiated from Python, but not very useful by\n\
itself.\n\
\n\
You can instantiate objects of this class by passing a set of construction\n\
parameters or another object of which the base type is ``FrameExtractor``.\n\
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
other\n\
  [FrameExtractor] an object of which is or inherits from a FrameExtractor\n\
  that will be deep-copied into a new instance.\n\
\n\
"
);

/**
 * Represents either an FrameExtractor
 */
typedef struct {
  PyObject_HEAD
  bob::ap::FrameExtractor* cxx;
} PyBobApFrameExtractorObject;

extern PyTypeObject PyBobApFrameExtractor_Type; //forward declaration

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
  static const char* const_kwlist[] = {"sampling_frequency, ", 0};
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

static int PyBobApFrameExtractor_InitParameters(PyBobApFrameExtractorObject* self, 
    PyObject *args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {
    "sampling_frequency", 
    "win_length_ms",
    "win_shift_ms",
    0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  double sampling_frequency = 0.;
  double win_length_ms = 20.;
  double win_shift_ms = 10.;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|dd", kwlist, &length)) return -1;

  try {
    self->cxx = new bob::ap::FrameExtractor(sampling_frequency, win_length_ms, win_shift_ms);
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
  return
# if PY_VERSION_HEX >= 0x03000000
  PyUnicode_FromFormat
# else
  PyString_FromFormat
# endif
  ("%s(sampling_frequency=%f, win_length_ms=%f, win_shift_ms=%f)", Py_TYPE(self)->tp_name, self->cxx->getSamplingFrequency(), self->cxx->getWinLengthMs(), self->cxx->getWinShiftMs());
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

  if (!PyNumber_Check(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' sampling frequency can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsFloat(o);
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

static PyGetSetDef PyBobApFrameExtractor_getseters[] = {
    {
      s_sampling_frequency_str,
      (getter)PyBobApFrameExtractor_GetSamplingFrequency,
      (setter)PyBobApFrameExtractor_SetSamplingFrequency,
      s_sampling_frequency_doc,
      0
    },
    {0}  /* Sentinel */
};

PyDoc_STRVAR(s_shape_str, "get_shape");
PyDoc_STRVAR(s_shape_doc,
"Computes the shape of the output features, given the size of an input\n\
array or an input array.\n\
");

static PyObject* PyBobApFrameExtractor_GetShape
(PyBobApFrameExtractorObject* self, void* /*closure*/) {
  return Py_BuildValue("(n)", self->cxx->getLength());
}

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
    (ternaryfunc)PyBobApFrameExtractor_Call,  /* tp_call */
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
    0,                                        /* tp_methods */
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
