/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @author Pavel Korshunov <Pavel.Korshunov@idiap.ch>
 * @date Thu  6 Feb 09:00:05 2014
 *
 * @brief Bindings to the base class bob::ap::Ceps
 */

#include <bob.blitz/cppapi.h>
#include <bob.blitz/cleanup.h>
#include <bob.extension/defines.h>
#include "types.h"

PyDoc_STRVAR(s_ceps_str, BOB_EXT_MODULE_PREFIX ".Ceps");

PyDoc_STRVAR(s_ceps_doc,
"Ceps(sampling_frequency, [win_length_ms=20., [win_shift_ms=10., [n_filters=24, [n_ceps=19, [f_min=0., [f_max=4000., [delta_win=2, [pre_emphasis_coeff=0.95, [mel_scale=True, [dct_norm=True, [normalize_mean=True, [rect_filter=False, [inverse_filter=False, [normalize_spectrum=False, [ssfc_features=False, [scfc_features=False, [scmc_features=False]]]]]]]]]]]]]]]]]) -> new Ceps\n\
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
normalize_mean\n\
  [bool] Tells whether frame should be normalized \n\
  by subtracting mean (True) or dividing by max_range (False)\n\
  ``True`` is the default value.\n\
\n\
rect_filter\n\
  [bool] tells whether to apply the filter in the\n\
  inversed order, i.e., from high frequencies to low\n\
  (set it to ``True''). ``False`` is the default value.\n\
\n\
inverse_filter\n\
  [bool] tells whether cepstral features are extracted\n\
  using a rectungular filter (set it to ``True``), i.e., RFCC features,\n\
  instead of the default filter (the default value is ``False``)\n\
\n\
normalize_spectrum\n\
  [bool] Tells whether to normalize the power spectrum of the signal.\n\
  The default value is ``False``.\n\
\n\
ssfc_features\n\
  [bool] Set to true if you want to compute\n\
  Subband Spectral Flux Coefficients (SSFC), which measures\n\
  the frame-by-frame change in the power spectrum\n\
\n\
scfc_features\n\
  [bool] Set to true if you want to compute\n\
  Spectral Centroid Frequency Coefficients (SCFC), which\n\
  capture detailed information about subbands similar to formant frequencies\n\
\n\
scmc_features\n\
  [bool] Set to true if you want to compute\n\
  Spectral Centroid Magnitude  Coefficients (SCMC), which\n\
  capture detailed information about subbands similar to SCFC features\n\
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

  o->parent.parent.parent.cxx = 0; // FrameExtractor
  o->parent.parent.cxx = 0; // Energy
  o->parent.cxx = 0; // Spectrogram
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
    "n_ceps",
    "f_min",
    "f_max",
    "delta_win",
    "pre_emphasis_coeff",
    "mel_scale",
    "dct_norm",
    "normalize_mean",
    "rect_filter",
    "inverse_filter",
    "normalize_spectrum",
    "ssfc_features",
    "scfc_features",
    "scmc_features",
    0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  double sampling_frequency = 0.;
  double win_length_ms = 20.;
  double win_shift_ms = 10.;
  Py_ssize_t n_filters = 24;
  Py_ssize_t n_ceps = 19;
  double f_min = 0.;
  double f_max = 8000.;
  Py_ssize_t delta_win = 2;
  double pre_emphasis_coeff = 0.95;
  PyObject* mel_scale = Py_True;
  PyObject* dct_norm = Py_False;
  PyObject* normalize_mean = Py_True;
  PyObject* rect_filter = Py_False;
  PyObject* inverse_filter = Py_False;
  PyObject* normalize_spectrum = Py_False;
  PyObject* ssfc_features = Py_False;
  PyObject* scfc_features = Py_False;
  PyObject* scmc_features = Py_False;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|ddnnddndOOOOOOOOO", kwlist,
        &sampling_frequency, &win_length_ms, &win_shift_ms, &n_filters,
        &n_ceps, &f_min, &f_max, &delta_win, &pre_emphasis_coeff,
        &mel_scale, &dct_norm, &normalize_mean, &rect_filter,
        &inverse_filter, &normalize_spectrum,
        &ssfc_features, &scfc_features, &scmc_features))
    return -1;

  bool mel_scale_ = PyObject_IsTrue(mel_scale);
  bool dct_norm_ = PyObject_IsTrue(dct_norm);
  bool normalize_mean_ = PyObject_IsTrue(normalize_mean);
  bool rect_filter_ = PyObject_IsTrue(rect_filter);
  bool inverse_filter_ = PyObject_IsTrue(inverse_filter);
  bool normalize_spectrum_ = PyObject_IsTrue(normalize_spectrum);
  bool ssfc_features_ = PyObject_IsTrue(ssfc_features);
  bool scfc_features_ = PyObject_IsTrue(scfc_features);
  bool scmc_features_ = PyObject_IsTrue(scmc_features);

  try {
    self->cxx = new bob::ap::Ceps(sampling_frequency,
        win_length_ms, win_shift_ms, n_filters, n_ceps, f_min, f_max,
        delta_win, pre_emphasis_coeff, mel_scale_, dct_norm_, normalize_mean_,
        rect_filter_, inverse_filter_, normalize_spectrum_,
        ssfc_features_, scfc_features_, scmc_features_);
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
  static const int MAXSIZE = 256;
  char buffer[MAXSIZE];
  Py_ssize_t n_filters = self->cxx->getNFilters();
  Py_ssize_t n_ceps = self->cxx->getNCeps();
  Py_ssize_t delta_win = self->cxx->getDeltaWin();
  auto count = std::snprintf(buffer, MAXSIZE, "%s(sampling_frequency=%f, win_length_ms=%f, win_shift_ms=%f, n_filters=%" PY_FORMAT_SIZE_T "d, n_ceps=%" PY_FORMAT_SIZE_T "d, f_min=%f, f_max=%f, delta_win=%" PY_FORMAT_SIZE_T "d, pre_emphasis_coeff=%f, mel_scale=%s, dct_norm=%s, normalize_mean=%s, rect_filter=%s, inverse_filter=%s, normalize_spectrum=%s, ssfc_features=%s, scfc_features=%s, scmc_features=%s)", Py_TYPE(self)->tp_name, self->cxx->getSamplingFrequency(), self->cxx->getWinLengthMs(), self->cxx->getWinShiftMs(), n_filters, n_ceps, self->cxx->getFMin(), self->cxx->getFMax(), delta_win, self->cxx->getPreEmphasisCoeff(), self->cxx->getMelScale()?"True":"False", self->cxx->getDctNorm()?"True":"False", self->cxx->getNormalizeMean()?"True":"False", self->cxx->getRectangularFilter()?"True":"False", self->cxx->getInverseFilter()?"True":"False", self->cxx->getNormalizeSpectrum()?"True":"False", self->cxx->getSSFCFeatures()?"True":"False", self->cxx->getSCFCFeatures()?"True":"False", self->cxx->getSCMCFeatures()?"True":"False");
  return
# if PY_VERSION_HEX >= 0x03000000
  PyUnicode_FromStringAndSize
# else
  PyString_FromStringAndSize
# endif
    (buffer, (count<=MAXSIZE)?count:MAXSIZE);
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

PyDoc_STRVAR(s_n_ceps_str, "n_ceps");
PyDoc_STRVAR(s_n_ceps_doc,
"The number of cepstral coefficients"
);

static PyObject* PyBobApCeps_GetNCeps
(PyBobApCepsObject* self, void* /*closure*/) {
  return Py_BuildValue("n", self->cxx->getNCeps());
}

static int PyBobApCeps_SetNCeps
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' n_ceps can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  Py_ssize_t n = PyNumber_AsSsize_t(o, PyExc_OverflowError);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setNCeps(n);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `n_ceps' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_delta_win_str, "delta_win");
PyDoc_STRVAR(s_delta_win_doc,
"The integer delta value used for computing the first and\n\
second order derivatives"
);

static PyObject* PyBobApCeps_GetDeltaWin
(PyBobApCepsObject* self, void* /*closure*/) {
  return Py_BuildValue("n", self->cxx->getDeltaWin());
}

static int PyBobApCeps_SetDeltaWin
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
    PyErr_Format(PyExc_TypeError, "`%s' delta_win can only be set using a number, not `%s'", Py_TYPE(self)->tp_name, Py_TYPE(o)->tp_name);
    return -1;
  }

  Py_ssize_t n = PyNumber_AsSsize_t(o, PyExc_OverflowError);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setDeltaWin(n);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `delta_win' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_dct_norm_str, "dct_norm");
PyDoc_STRVAR(s_dct_norm_doc,
"A factor by which the cepstral coefficients are multiplied"
);

static PyObject* PyBobApCeps_GetDctNorm
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getDctNorm()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetDctNorm
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setDctNorm(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `dct_norm' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_with_energy_str, "with_energy");
PyDoc_STRVAR(s_with_energy_doc,
"Tells if we add the energy to the output feature"
);

static PyObject* PyBobApCeps_GetWithEnergy
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getWithEnergy()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetWithEnergy
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setWithEnergy(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `with_energy' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_with_delta_str, "with_delta");
PyDoc_STRVAR(s_with_delta_doc,
"Tells if we add the first derivatives to the output feature"
);

static PyObject* PyBobApCeps_GetWithDelta
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getWithDelta()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetWithDelta
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setWithDelta(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `with_delta' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_with_delta_delta_str, "with_delta_delta");
PyDoc_STRVAR(s_with_delta_delta_doc,
"Tells if we add the second derivatives to the output feature"
);

static PyObject* PyBobApCeps_GetWithDeltaDelta
(PyBobApCepsObject* self, void* /*closure*/) {
  if (self->cxx->getWithDeltaDelta()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApCeps_SetWithDeltaDelta
(PyBobApCepsObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setWithDeltaDelta(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `with_delta_delta' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

static PyGetSetDef PyBobApCeps_getseters[] = {
    {
      s_n_ceps_str,
      (getter)PyBobApCeps_GetNCeps,
      (setter)PyBobApCeps_SetNCeps,
      s_n_ceps_doc,
      0
    },
    {
      s_delta_win_str,
      (getter)PyBobApCeps_GetDeltaWin,
      (setter)PyBobApCeps_SetDeltaWin,
      s_delta_win_doc,
      0
    },
    {
      s_dct_norm_str,
      (getter)PyBobApCeps_GetDctNorm,
      (setter)PyBobApCeps_SetDctNorm,
      s_dct_norm_doc,
      0
    },
    {
      s_with_energy_str,
      (getter)PyBobApCeps_GetWithEnergy,
      (setter)PyBobApCeps_SetWithEnergy,
      s_with_energy_doc,
      0
    },
    {
      s_with_delta_str,
      (getter)PyBobApCeps_GetWithDelta,
      (setter)PyBobApCeps_SetWithDelta,
      s_with_delta_doc,
      0
    },
    {
      s_with_delta_delta_str,
      (getter)PyBobApCeps_GetWithDeltaDelta,
      (setter)PyBobApCeps_SetWithDeltaDelta,
      s_with_delta_delta_doc,
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

  return PyBlitzArray_NUMPY_WRAP(Py_BuildValue("O", output));

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
