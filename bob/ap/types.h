/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Thu 13 Feb 2014 15:07:49 CET
 *
 * @brief All types that need cross-compilation on different units
 */

#ifndef BOB_AP_TYPES_H
#define BOB_AP_TYPES_H

#include <bob.ap/FrameExtractor.h>
#include <bob.ap/Energy.h>
#include <bob.ap/Spectrogram.h>
#include <bob.ap/Ceps.h>

/**
 * Represents either an FrameExtractor
 */
typedef struct {
  PyObject_HEAD
  bob::ap::FrameExtractor* cxx;
} PyBobApFrameExtractorObject;

extern PyTypeObject PyBobApFrameExtractor_Type; //forward declaration

/**
 * Represents either the Energy extractor
 */
typedef struct {
  PyBobApFrameExtractorObject parent;
  bob::ap::Energy* cxx;
} PyBobApEnergyObject;

extern PyTypeObject PyBobApEnergy_Type; //forward declaration

/**
 * Represents either the Spectrogram extractor
 */
typedef struct {
  PyBobApEnergyObject parent;
  bob::ap::Spectrogram* cxx;
} PyBobApSpectrogramObject;

extern PyTypeObject PyBobApSpectrogram_Type; //forward declaration

/**
 * Represents either the Ceps extractor
 */
typedef struct {
  PyBobApSpectrogramObject parent;
  bob::ap::Ceps* cxx;
} PyBobApCepsObject;

extern PyTypeObject PyBobApCeps_Type; //forward declaration

#endif /* BOB_AP_TYPES_H */
