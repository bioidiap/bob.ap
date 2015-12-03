/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Fri 25 Oct 16:54:55 2013
 *
 * @brief Bindings to bob::ap
 */

#ifdef NO_IMPORT_ARRAY
#undef NO_IMPORT_ARRAY
#endif
#include <bob.blitz/capi.h>
#include <bob.blitz/cleanup.h>
#include <bob.core/api.h>
#include <bob.sp/api.h>

#include "types.h"

static PyMethodDef module_methods[] = {
    {0}  /* Sentinel */
};

PyDoc_STRVAR(module_docstr, "Bob audio processing toolkit");

#if PY_VERSION_HEX >= 0x03000000
static PyModuleDef module_definition = {
  PyModuleDef_HEAD_INIT,
  BOB_EXT_MODULE_NAME,
  module_docstr,
  -1,
  module_methods,
  0, 0, 0, 0
};
#endif

static PyObject* create_module (void) {

  PyBobApFrameExtractor_Type.tp_new = PyType_GenericNew;
  if (PyType_Ready(&PyBobApFrameExtractor_Type) < 0) return 0;

  PyBobApEnergy_Type.tp_base = &PyBobApFrameExtractor_Type;
  if (PyType_Ready(&PyBobApEnergy_Type) < 0) return 0;

  PyBobApSpectrogram_Type.tp_base = &PyBobApEnergy_Type;
  if (PyType_Ready(&PyBobApSpectrogram_Type) < 0) return 0;

  PyBobApCeps_Type.tp_base = &PyBobApSpectrogram_Type;
  if (PyType_Ready(&PyBobApCeps_Type) < 0) return 0;

# if PY_VERSION_HEX >= 0x03000000
  PyObject* m = PyModule_Create(&module_definition);
  auto m_ = make_xsafe(m);
  const char* ret = "O";
# else
  PyObject* m = Py_InitModule3(BOB_EXT_MODULE_NAME, module_methods, module_docstr);
  const char* ret = "N";
# endif
  if (!m) return 0;

  /* register the types to python */
  Py_INCREF(&PyBobApFrameExtractor_Type);
  if (PyModule_AddObject(m, "FrameExtractor", (PyObject *)&PyBobApFrameExtractor_Type) < 0) return 0;

  Py_INCREF(&PyBobApEnergy_Type);
  if (PyModule_AddObject(m, "Energy", (PyObject *)&PyBobApEnergy_Type) < 0) return 0;

  Py_INCREF(&PyBobApSpectrogram_Type);
  if (PyModule_AddObject(m, "Spectrogram", (PyObject *)&PyBobApSpectrogram_Type) < 0) return 0;

  Py_INCREF(&PyBobApCeps_Type);
  if (PyModule_AddObject(m, "Ceps", (PyObject *)&PyBobApCeps_Type) < 0) return 0;

  /* imports dependencies */
  if (import_bob_blitz() < 0) return 0;
  if (import_bob_core_logging() < 0) return 0;
  if (import_bob_sp() < 0) return 0;

  return Py_BuildValue(ret, m);
}

PyMODINIT_FUNC BOB_EXT_ENTRY_NAME (void) {
# if PY_VERSION_HEX >= 0x03000000
  return
# endif
    create_module();
}
