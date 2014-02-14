/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Thu  6 Feb 09:00:05 2014
 *
 * @brief Bindings to the base class bob::ap::Ceps
 */

#include <xbob.blitz/cppapi.h>
#include <xbob.blitz/cleanup.h>
#include "types.h"

PyDoc_STRVAR(s_ceps_str, XBOB_EXT_MODULE_PREFIX ".Ceps");

PyDoc_STRVAR(s_ceps_doc,
"Ceps(sampling_frequency, [win_length_ms=20., [win_shift_ms=10., [n_filters=24, [n_ceps=19, [f_min=0., [f_max=4000., [delta_win=2, [pre_emphasis_coeff=0.95, [mel_scale=True, [dct_norm=True]]]]]]]]]]) -> new Ceps\n\
Ceps(other) -> new Ceps\n\
\n\
Objects of this class, after configuration, can extract the\n\
cepstral coefficients from 1D audio array/signals.\n\
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
n_filters\n\
  [int] the number of filter bands\n\
\n\
n_ceps\n\
  [int] the number of cepstral coefficients\n\
\n\
f_min\n\
  [double] the minimum frequency of the filter bank\n\
\n\
f_max\n\
  [double] the maximum frequency of the filter bank\n\
\n\
delta_win\n\
  [int] The integer delta value used for computing the\n\
  first and second order derivatives\n\
\n\
pre_emphasis_coeff\n\
  [double] the coefficient used for the pre-emphasis\n\
\n\
mel_scale\n\
  [bool] tells whether cepstral features are extracted\n\
  on a linear (LFCC, set it to ``False``) or Mel (MFCC,\n\
  set it to ``True`` - the default)\n\
\n\
dct_norm\n\
  [bool] A factor by which the cepstral coefficients are\n\
  multiplied\n\
\n\
\n\
other\n\
  [Ceps] an object of which is or inherits from ``Ceps``\n\
  that will be deep-copied into a new instance.\n\
\n\
"
);

int PyBobApCeps_Check(PyObject* o) {
  return PyObject_IsInstance(o, reinterpret_cast<PyObject*>(&PyBobApCeps_Type));
}

static void PyBobApCeps_Delete (PyBobApCepsObject* o) {

  o->parent.parent.parent.cxx = 0;
  o->parent.parent.cxx = 0;
  o->parent.cxx = 0;
  delete o->cxx;
  Py_TYPE(o)->tp_free((PyObject*)o);

}

static int PyBobApCeps_InitCopy
(PyBobApCepsObject* self, PyObject* args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"other", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  PyObject* other = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
        &PyBobApCeps_Type, &other)) return -1;

  auto copy = reinterpret_cast<PyBobApCepsObject*>(other);

  try {
    self->cxx = new bob::ap::Ceps(*(copy->cxx));
    if (!self->cxx) {
      PyErr_Format(PyExc_MemoryError, "cannot create new object of type `%s' - no more memory", Py_TYPE(self)->tp_name);
      return -1;
    }
    self->parent.parent.parent.cxx = self->cxx;
    self->parent.parent.cxx = self->cxx;
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

static int PyBobApCeps_InitParameters
(PyBobApCepsObject* self, PyObject *args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {
    "sampling_frequency",
    "win_length_ms",
    "win_shift_ms",
    "n_filters",
    "f_min",
    "f_max",
    "pre_emphasis_coeff",
    "mel_scale",
    0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  double sampling_frequency = 0.;
  double win_length_ms = 20.;
  double win_shift_ms = 10.;
  Py_ssize_t n_filters = 24;
  double f_min = 0.;
  double f_max = 0.;
  double pre_emphasis_coeff = 0.;
  PyObject* mel_scale = Py_True;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|ddndddO", kwlist,
        &sampling_frequency, &win_length_ms, &win_shift_ms,
        &n_filters, &f_min, &f_max, &pre_emphasis_coeff, &mel_scale))
    return -1;

  bool mel_scale_ = PyObject_IsTrue(mel_scale);

  try {
    self->cxx = new bob::ap::Ceps(sampling_frequency,
        win_length_ms, win_shift_ms, n_filters, f_min, f_max,
        pre_emphasis_coeff, mel_scale_);
    if (!self->cxx) {
      PyErr_Format(PyExc_MemoryError, "cannot create new object of type `%s' - no more memory", Py_TYPE(self)->tp_name);
      return -1;
    }
    self->parent.parent.cxx = self->cxx;
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

static int PyBobApCeps_Init(PyBobApCepsObject* self,
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

        if (PyBobApCeps_Check(arg)) {
          return PyBobApCeps_InitCopy(self, args, kwds);
        }

        else {
          return PyBobApCeps_InitParameters(self, args, kwds);
        }

        PyErr_Format(PyExc_TypeError, "cannot initialize `%s' with `%s' (see help)", Py_TYPE(self)->tp_name, Py_TYPE(arg)->tp_name);

      }

      break;

    default:

      return PyBobApCeps_InitParameters(self, args, kwds);

  }

  return -1;

}

static PyObject* PyBobApCeps_Repr(PyBobApCepsObject* self) {
  Py_ssize_t n_filters = self->cxx->getNFilters();
  return
# if PY_VERSION_HEX >= 0x03000000
  PyUnicode_FromFormat
# else
  PyString_FromFormat
# endif
  ("%s(sampling_frequency=%f, win_length_ms=%f, win_shift_ms=%f, n_filters=%" PY_FORMAT_SIZE_T "d, f_min=%f, f_max=%f, pre_emphasis_coeff=%f, mel_scale=%s)", Py_TYPE(self)->tp_name, self->cxx->getSamplingFrequency(), self->cxx->getWinLengthMs(), self->cxx->getWinShiftMs(), n_filters, self->cxx->getFMin(), self->cxx->getFMax(), self->cxx->getPreEmphasisCoeff(), self->cxx->getMelScale()?"True":"False");
}

static PyObject* PyBobApCeps_RichCompare (PyBobApCepsObject* self,
    PyObject* other, int op) {

  if (!PyBobApCeps_Check(other)) {
    PyErr_Format(PyExc_TypeError, "cannot compare `%s' with `%s'",
        Py_TYPE(self)->tp_name, Py_TYPE(other)->tp_name);
    return 0;
  }

  auto other_ = reinterpret_cast<PyBobApCepsObject*>(other);

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

PyDoc_STRVAR(s_n_filters_str, "n_filters");
PyDoc_STRVAR(s_n_filters_doc,
"The number of filter bands"
);

static PyObject* PyBobApCeps_GetNFilters
(PyBobApCepsObject* self, void* /*closure*/) {
  return Py_BuildValue("n", self->cxx->getNFilters());
}

static int PyBobApCeps_SetNFilters
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  if (!PyNumber_Check(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' n_filters can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  Py_ssize_t n = PyNumber_AsSsize_t(o, PyExc_OverflowError);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setNFilters(n);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `n_filters' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_f_min_str, "f_min");
PyDoc_STRVAR(s_f_min_doc,
"The minimum frequency of the filter bank"
);

static PyObject* PyBobApCeps_GetFMin
(PyBobApCepsObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getFMin());
}

static int PyBobApCeps_SetFMin
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  if (!PyNumber_Check(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' f_min can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setFMin(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `f_min' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_f_max_str, "f_max");
PyDoc_STRVAR(s_f_max_doc,
"The maximum frequency of the filter bank"
);

static PyObject* PyBobApCeps_GetFMax
(PyBobApCepsObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getFMax());
}

static int PyBobApCeps_SetFMax
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  if (!PyNumber_Check(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' f_max can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setFMax(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `f_max' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_pre_emphasis_coeff_str, "pre_emphasis_coeff");
PyDoc_STRVAR(s_pre_emphasis_coeff_doc,
"The coefficient used for the pre-emphasis"
);

static PyObject* PyBobApCeps_GetPreEmphasisCoeff
(PyBobApCepsObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getPreEmphasisCoeff());
}

static int PyBobApCeps_SetPreEmphasisCoeff
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  if (!PyNumber_Check(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' pre_emphasis_coeff can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  double d = PyFloat_AsDouble(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setPreEmphasisCoeff(d);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `pre_emphasis_coeff' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_mel_scale_str, "mel_scale");
PyDoc_STRVAR(s_mel_scale_doc,
"Tells whether cepstral features are extracted on a linear\n\
(LFCC) or Mel (MFCC) scale\n\
");

static PyObject* PyBobApCeps_GetMelScale
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getMelScale()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetMelScale
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setMelScale(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `mel_scale' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_energy_filter_str, "energy_filter");
PyDoc_STRVAR(s_energy_filter_doc,
"Tells whether we use the energy or the square root of the energy"
);

static PyObject* PyBobApCeps_GetEnergyFilter
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getEnergyFilter()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetEnergyFilter
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setEnergyFilter(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `energy_filter' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_log_filter_str, "log_filter");
PyDoc_STRVAR(s_log_filter_doc,
"Tells whether we use the log triangular filter or the triangular filter"
);

static PyObject* PyBobApCeps_GetLogFilter
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getLogFilter()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetLogFilter
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setLogFilter(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `log_filter' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_energy_bands_str, "energy_bands");
PyDoc_STRVAR(s_energy_bands_doc,
"Tells whether we compute a ceps or energy bands"
);

static PyObject* PyBobApCeps_GetEnergyBands
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getEnergyBands()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetEnergyBands
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setEnergyBands(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `energy_bands' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

static PyGetSetDef PyBobApCeps_getseters[] = {
    {
      s_n_filters_str,
      (getter)PyBobApCeps_GetNFilters,
      (setter)PyBobApCeps_SetNFilters,
      s_n_filters_doc,
      0
    },
    {
      s_f_min_str,
      (getter)PyBobApCeps_GetFMin,
      (setter)PyBobApCeps_SetFMin,
      s_f_min_doc,
      0
    },
    {
      s_f_max_str,
      (getter)PyBobApCeps_GetFMax,
      (setter)PyBobApCeps_SetFMax,
      s_f_max_doc,
      0
    },
    {
      s_pre_emphasis_coeff_str,
      (getter)PyBobApCeps_GetPreEmphasisCoeff,
      (setter)PyBobApCeps_SetPreEmphasisCoeff,
      s_pre_emphasis_coeff_doc,
      0
    },
    {
      s_mel_scale_str,
      (getter)PyBobApCeps_GetMelScale,
      (setter)PyBobApCeps_SetMelScale,
      s_mel_scale_doc,
      0
    },
    {
      s_energy_filter_str,
      (getter)PyBobApCeps_GetEnergyFilter,
      (setter)PyBobApCeps_SetEnergyFilter,
      s_energy_filter_doc,
      0
    },
    {
      s_log_filter_str,
      (getter)PyBobApCeps_GetLogFilter,
      (setter)PyBobApCeps_SetLogFilter,
      s_log_filter_doc,
      0
    },
    {
      s_energy_bands_str,
      (getter)PyBobApCeps_GetEnergyBands,
      (setter)PyBobApCeps_SetEnergyBands,
      s_energy_bands_doc,
      0
    },
    {0}  /* Sentinel */
};

static PyObject* PyBobApCeps_Call
(PyBobApCepsObject* self, PyObject *args, PyObject* kwds) {

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
      PyErr_Format(PyExc_TypeError, "`%s' only supports 2D 64-bit float arrays for output array `output'", Py_TYPE(self)->tp_name);
      return 0;
    }

    if (output->ndim != 2) {
      PyErr_Format(PyExc_TypeError, "`%s' only supports 2D 64-bit float arrays for output array `output'", Py_TYPE(self)->tp_name);
      return 0;
    }

  }

  else {

    Py_ssize_t length[2];
    auto s = self->cxx->getShape(*bz_input);
    length[0] = s(0);
    length[1] = s(1);
    output = (PyBlitzArrayObject*)PyBlitzArray_SimpleNew(NPY_FLOAT64, 2, length);
    if (!output) return 0;
    output_ = make_safe(output);

  }

  auto bz_output = PyBlitzArrayCxx_AsBlitz<double,2>(output);

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

  Py_INCREF(output);
  return PyBlitzArray_NUMPY_WRAP(reinterpret_cast<PyObject*>(output));

}

PyTypeObject PyBobApCeps_Type = {
    PyVarObject_HEAD_INIT(0, 0)
    s_ceps_str,                               /*tp_name*/
    sizeof(PyBobApCepsObject),                /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)PyBobApCeps_Delete,           /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    (reprfunc)PyBobApCeps_Repr,               /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash */
    (ternaryfunc)PyBobApCeps_Call,            /* tp_call */
    (reprfunc)PyBobApCeps_Repr,               /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    s_ceps_doc,                               /* tp_doc */
    0,		                                    /* tp_traverse */
    0,		                                    /* tp_clear */
    (richcmpfunc)PyBobApCeps_RichCompare,     /* tp_richcompare */
    0,		                                    /* tp_weaklistoffset */
    0,		                                    /* tp_iter */
    0,		                                    /* tp_iternext */
    0,                                        /* tp_methods */
    0,                                        /* tp_members */
    PyBobApCeps_getseters,                    /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    (initproc)PyBobApCeps_Init,               /* tp_init */
    0,                                        /* tp_alloc */
    0,                                        /* tp_new */
};
