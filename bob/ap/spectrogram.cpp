/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @author Pavel Korshunov <Pavel.Korshunov@idiap.ch>
 * @date Thu  6 Feb 09:00:05 2014
 *
 * @brief Bindings to the base class bob::ap::Spectrogram
 */

#include <bob.blitz/cppapi.h>
#include <bob.blitz/cleanup.h>
#include <bob.extension/defines.h>
#include "types.h"

PyDoc_STRVAR(s_spectrogram_str, BOB_EXT_MODULE_PREFIX ".Spectrogram");

PyDoc_STRVAR(s_spectrogram_doc,
"Spectrogram(sampling_frequency, [win_length_ms=20., [win_shift_ms=10., [n_filters=24, [f_min=0., [f_max=4000., [pre_emphasis_coeff=0.95, [mel_scale=True, [normalize_mean=True, [rect_filter=False, [inverse_filter=False, [normalize_spectrum=False, [ssfc_features=False, [scfc_features=False, [scmc_features=False]]]]]]]]]]]]]]) -> new Spectrogram\n\
Spectrogram(other) -> new Spectrogram\n\
\n\
Objects of this class, after configuration, can extract the\n\
spectrogram from 1D audio array/signals.\n\
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
f_min\n\
  [double] the minimum frequency of the filter bank\n\
\n\
f_max\n\
  [double] the maximum frequency of the filter bank\n\
\n\
pre_emphasis_coeff\n\
  [double] the coefficient used for the pre-emphasis\n\
\n\
mel_scale\n\
  [bool] tells whether cepstral features are extracted\n\
  on a linear (LFCC, set it to ``False``) or Mel (MFCC,\n\
  set it to ``True`` - the default)\n\
\n\
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
  [Spectrogram] an object of which is or inherits from ``Spectrogram``\n\
  that will be deep-copied into a new instance.\n\
\n\
"
);

int PyBobApSpectrogram_Check(PyObject* o) {
  return PyObject_IsInstance(o, reinterpret_cast<PyObject*>(&PyBobApSpectrogram_Type));
}

static void PyBobApSpectrogram_Delete (PyBobApSpectrogramObject* o) {

  o->parent.parent.cxx = 0;
  o->parent.cxx = 0;
  delete o->cxx;
  Py_TYPE(o)->tp_free((PyObject*)o);

}

static int PyBobApSpectrogram_InitCopy
(PyBobApSpectrogramObject* self, PyObject* args, PyObject* kwds) {

  /* Parses input arguments in a single shot */
  static const char* const_kwlist[] = {"other", 0};
  static char** kwlist = const_cast<char**>(const_kwlist);

  PyObject* other = 0;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist,
        &PyBobApSpectrogram_Type, &other)) return -1;

  auto copy = reinterpret_cast<PyBobApSpectrogramObject*>(other);

  try {
    self->cxx = new bob::ap::Spectrogram(*(copy->cxx));
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

  return 0;

}

static int PyBobApSpectrogram_InitParameters
(PyBobApSpectrogramObject* self, PyObject *args, PyObject* kwds) {

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
  double f_min = 0.;
  double f_max = 8000.;
  double pre_emphasis_coeff = 0.95;
  PyObject* mel_scale = Py_True;
  PyObject* normalize_mean = Py_True;
  PyObject* rect_filter = Py_False;
  PyObject* inverse_filter = Py_False;
  PyObject* normalize_spectrum = Py_False;
  PyObject* ssfc_features = Py_False;
  PyObject* scfc_features = Py_False;
  PyObject* scmc_features = Py_False;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|ddndddOOOOOOOO", kwlist,
        &sampling_frequency, &win_length_ms, &win_shift_ms,
        &n_filters, &f_min, &f_max, &pre_emphasis_coeff, &mel_scale,
        &normalize_mean, &rect_filter, &inverse_filter, &normalize_spectrum,
        &ssfc_features, &scfc_features, &scmc_features))
    return -1;

  bool mel_scale_ = PyObject_IsTrue(mel_scale);
  bool normalize_mean_ = PyObject_IsTrue(normalize_mean);
  bool rect_filter_ = PyObject_IsTrue(rect_filter);
  bool inverse_filter_ = PyObject_IsTrue(inverse_filter);
  bool normalize_spectrum_ = PyObject_IsTrue(normalize_spectrum);
  bool ssfc_features_ = PyObject_IsTrue(ssfc_features);
  bool scfc_features_ = PyObject_IsTrue(scfc_features);
  bool scmc_features_ = PyObject_IsTrue(scmc_features);

  try {
    self->cxx = new bob::ap::Spectrogram(sampling_frequency,
        win_length_ms, win_shift_ms, n_filters, f_min, f_max,
        pre_emphasis_coeff, mel_scale_, normalize_mean_, rect_filter_,
        inverse_filter_, normalize_spectrum_,
        ssfc_features_, scfc_features_, scmc_features_);
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

static int PyBobApSpectrogram_Init(PyBobApSpectrogramObject* self,
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

        if (PyBobApSpectrogram_Check(arg)) {
          return PyBobApSpectrogram_InitCopy(self, args, kwds);
        }

        else {
          return PyBobApSpectrogram_InitParameters(self, args, kwds);
        }

        PyErr_Format(PyExc_TypeError, "cannot initialize `%s' with `%s' (see help)", Py_TYPE(self)->tp_name, Py_TYPE(arg)->tp_name);

      }

      break;

    default:

      return PyBobApSpectrogram_InitParameters(self, args, kwds);

  }

  return -1;

}

static PyObject* PyBobApSpectrogram_Repr(PyBobApSpectrogramObject* self) {
  static const int MAXSIZE = 256;
  char buffer[MAXSIZE];
  Py_ssize_t n_filters = self->cxx->getNFilters();
  auto count = std::snprintf(buffer, MAXSIZE, "%s(sampling_frequency=%f, win_length_ms=%f, win_shift_ms=%f,n_filters=%" PY_FORMAT_SIZE_T "d, f_min=%f, f_max=%f, pre_emphasis_coeff=%f, mel_scale=%s,  normalize_mean=%s, rect_filter=%s, inverse_filter=%s, normalize_spectrum=%s, ssfc_features=%s, scfc_features=%s, scmc_features=%s)", Py_TYPE(self)->tp_name, self->cxx->getSamplingFrequency(), self->cxx->getWinLengthMs(), self->cxx->getWinShiftMs(), n_filters, self->cxx->getFMin(), self->cxx->getFMax(), self->cxx->getPreEmphasisCoeff(), self->cxx->getMelScale()?"True":"False", self->cxx->getNormalizeMean()?"True":"False", self->cxx->getRectangularFilter()?"True":"False", self->cxx->getInverseFilter()?"True":"False", self->cxx->getNormalizeSpectrum()?"True":"False", self->cxx->getSSFCFeatures()?"True":"False", self->cxx->getSCFCFeatures()?"True":"False", self->cxx->getSCMCFeatures()?"True":"False");
  return
# if PY_VERSION_HEX >= 0x03000000
  PyUnicode_FromStringAndSize
# else
  PyString_FromStringAndSize
# endif
    (buffer, (count<=MAXSIZE)?count:MAXSIZE);
}

static PyObject* PyBobApSpectrogram_RichCompare (PyBobApSpectrogramObject* self,
    PyObject* other, int op) {

  if (!PyBobApSpectrogram_Check(other)) {
    PyErr_Format(PyExc_TypeError, "cannot compare `%s' with `%s'",
        Py_TYPE(self)->tp_name, Py_TYPE(other)->tp_name);
    return 0;
  }

  auto other_ = reinterpret_cast<PyBobApSpectrogramObject*>(other);

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

static PyObject* PyBobApSpectrogram_GetNFilters
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  return Py_BuildValue("n", self->cxx->getNFilters());
}

static int PyBobApSpectrogram_SetNFilters
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
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

static PyObject* PyBobApSpectrogram_GetFMin
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getFMin());
}

static int PyBobApSpectrogram_SetFMin
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
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

static PyObject* PyBobApSpectrogram_GetFMax
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getFMax());
}

static int PyBobApSpectrogram_SetFMax
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
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

static PyObject* PyBobApSpectrogram_GetPreEmphasisCoeff
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  return Py_BuildValue("d", self->cxx->getPreEmphasisCoeff());
}

static int PyBobApSpectrogram_SetPreEmphasisCoeff
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  if (!PyBob_NumberCheck(o)) {
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

static PyObject* PyBobApSpectrogram_GetMelScale
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getMelScale()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetMelScale
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

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


PyDoc_STRVAR(s_rect_filter_str, "rect_filter");
PyDoc_STRVAR(s_rect_filter_doc,
"Tells whether cepstral features are extracted using a rectangular scale\n\
");

static PyObject* PyBobApSpectrogram_GetRectangularFilter
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getRectangularFilter()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetRectangularFilter
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setRectangularFilter(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `rect_filter' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;
}


PyDoc_STRVAR(s_inverse_filter_str, "inverse_filter");
PyDoc_STRVAR(s_inverse_filter_doc,
"Tells whether the filter is applied in the inversed order when cepstral features are extracted\n\
");

static PyObject* PyBobApSpectrogram_GetInverseFilter
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getInverseFilter()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetInverseFilter
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setInverseFilter(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `inverse_filter' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;
}


PyDoc_STRVAR(s_normalize_spectrum_str, "normalize_spectrum");
PyDoc_STRVAR(s_normalize_spectrum_doc,
"Tells whether the filter is applied in the inversed order when cepstral features are extracted\n\
");

static PyObject* PyBobApSpectrogram_GetNormalizeSpectrum
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getNormalizeSpectrum()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetNormalizeSpectrum
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setNormalizeSpectrum(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `normalize_spectrum' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;
}


PyDoc_STRVAR(s_ssfc_features_str, "ssfc_features");
PyDoc_STRVAR(s_ssfc_features_doc,
"Make true if you want to compute SSFC features"
);

static PyObject* PyBobApSpectrogram_GetSSFCFeatures
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getSSFCFeatures()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetSSFCFeatures
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setSSFCFeatures(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `ssfc_features' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;

}

PyDoc_STRVAR(s_scfc_features_str, "scfc_features");
PyDoc_STRVAR(s_scfc_features_doc,
"Make true if you want to compute SCFC features"
);

static PyObject* PyBobApSpectrogram_GetSCFCFeatures
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getSCFCFeatures()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetSCFCFeatures
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setSCFCFeatures(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `scfc_features' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;
}

PyDoc_STRVAR(s_scmc_features_str, "scmc_features");
PyDoc_STRVAR(s_scmc_features_doc,
"Make true if you want to compute SCMC features"
);

static PyObject* PyBobApSpectrogram_GetSCMCFeatures
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getSCMCFeatures()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetSCMCFeatures
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

  bool b = PyObject_IsTrue(o);
  if (PyErr_Occurred()) return -1;

  try {
    self->cxx->setSCMCFeatures(b);
  }
  catch (std::exception& ex) {
    PyErr_SetString(PyExc_RuntimeError, ex.what());
    return -1;
  }
  catch (...) {
    PyErr_Format(PyExc_RuntimeError, "cannot reset `scmc_features' of %s: unknown exception caught", Py_TYPE(self)->tp_name);
    return -1;
  }

  return 0;
}


PyDoc_STRVAR(s_energy_filter_str, "energy_filter");
PyDoc_STRVAR(s_energy_filter_doc,
"Tells whether we use the energy or the square root of the energy"
);

static PyObject* PyBobApSpectrogram_GetEnergyFilter
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getEnergyFilter()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetEnergyFilter
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

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

static PyObject* PyBobApSpectrogram_GetLogFilter
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getLogFilter()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetLogFilter
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

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
"Tells whether we compute a spectrogram or energy bands"
);

static PyObject* PyBobApSpectrogram_GetEnergyBands
(PyBobApSpectrogramObject* self, void* /*closure*/) {
  if (self->cxx->getEnergyBands()) Py_RETURN_TRUE;
  else Py_RETURN_FALSE;
}

static int PyBobApSpectrogram_SetEnergyBands
(PyBobApSpectrogramObject* self, PyObject* o, void* /*closure*/) {

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

static PyGetSetDef PyBobApSpectrogram_getseters[] = {
    {
      s_n_filters_str,
      (getter)PyBobApSpectrogram_GetNFilters,
      (setter)PyBobApSpectrogram_SetNFilters,
      s_n_filters_doc,
      0
    },
    {
      s_f_min_str,
      (getter)PyBobApSpectrogram_GetFMin,
      (setter)PyBobApSpectrogram_SetFMin,
      s_f_min_doc,
      0
    },
    {
      s_f_max_str,
      (getter)PyBobApSpectrogram_GetFMax,
      (setter)PyBobApSpectrogram_SetFMax,
      s_f_max_doc,
      0
    },
    {
      s_pre_emphasis_coeff_str,
      (getter)PyBobApSpectrogram_GetPreEmphasisCoeff,
      (setter)PyBobApSpectrogram_SetPreEmphasisCoeff,
      s_pre_emphasis_coeff_doc,
      0
    },
    {
      s_mel_scale_str,
      (getter)PyBobApSpectrogram_GetMelScale,
      (setter)PyBobApSpectrogram_SetMelScale,
      s_mel_scale_doc,
      0
    },
    {
      s_rect_filter_str,
      (getter)PyBobApSpectrogram_GetRectangularFilter,
      (setter)PyBobApSpectrogram_SetRectangularFilter,
      s_rect_filter_doc,
      0
    },
    {
      s_inverse_filter_str,
      (getter)PyBobApSpectrogram_GetInverseFilter,
      (setter)PyBobApSpectrogram_SetInverseFilter,
      s_inverse_filter_doc,
      0
    },
    {
      s_normalize_spectrum_str,
      (getter)PyBobApSpectrogram_GetNormalizeSpectrum,
      (setter)PyBobApSpectrogram_SetNormalizeSpectrum,
      s_normalize_spectrum_doc,
      0
    },
    {
      s_ssfc_features_str,
      (getter)PyBobApSpectrogram_GetSSFCFeatures,
      (setter)PyBobApSpectrogram_SetSSFCFeatures,
      s_ssfc_features_doc,
      0
    },
    {
      s_scfc_features_str,
      (getter)PyBobApSpectrogram_GetSCFCFeatures,
      (setter)PyBobApSpectrogram_SetSCFCFeatures,
      s_scfc_features_doc,
      0
    },
    {
      s_scmc_features_str,
      (getter)PyBobApSpectrogram_GetSCMCFeatures,
      (setter)PyBobApSpectrogram_SetSCMCFeatures,
      s_scmc_features_doc,
      0
    },
    {
      s_energy_filter_str,
      (getter)PyBobApSpectrogram_GetEnergyFilter,
      (setter)PyBobApSpectrogram_SetEnergyFilter,
      s_energy_filter_doc,
      0
    },
    {
      s_log_filter_str,
      (getter)PyBobApSpectrogram_GetLogFilter,
      (setter)PyBobApSpectrogram_SetLogFilter,
      s_log_filter_doc,
      0
    },
    {
      s_energy_bands_str,
      (getter)PyBobApSpectrogram_GetEnergyBands,
      (setter)PyBobApSpectrogram_SetEnergyBands,
      s_energy_bands_doc,
      0
    },
    {0}  /* Sentinel */
};

static PyObject* PyBobApSpectrogram_Call
(PyBobApSpectrogramObject* self, PyObject *args, PyObject* kwds) {

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

PyTypeObject PyBobApSpectrogram_Type = {
    PyVarObject_HEAD_INIT(0, 0)
    s_spectrogram_str,                        /*tp_name*/
    sizeof(PyBobApSpectrogramObject),         /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)PyBobApSpectrogram_Delete,    /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    (reprfunc)PyBobApSpectrogram_Repr,        /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash */
    (ternaryfunc)PyBobApSpectrogram_Call,     /* tp_call */
    (reprfunc)PyBobApSpectrogram_Repr,        /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    s_spectrogram_doc,                        /* tp_doc */
    0,		                                    /* tp_traverse */
    0,		                                    /* tp_clear */
    (richcmpfunc)PyBobApSpectrogram_RichCompare,    /* tp_richcompare */
    0,		                                    /* tp_weaklistoffset */
    0,		                                    /* tp_iter */
    0,		                                    /* tp_iternext */
    0,                                        /* tp_methods */
    0,                                        /* tp_members */
    PyBobApSpectrogram_getseters,             /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    (initproc)PyBobApSpectrogram_Init,        /* tp_init */
    0,                                        /* tp_alloc */
    0,                                        /* tp_new */
};
