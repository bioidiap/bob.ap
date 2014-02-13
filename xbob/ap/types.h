/**
 * @author Andre Anjos <andre.anjos@idiap.ch>
 * @date Thu 13 Feb 2014 15:07:49 CET
 *
 * @brief All types that need cross-compilation on different units
 */

#ifndef XBOB_AP_TYPES_H
#define XBOB_AP_TYPES_H

#include <bob/ap/FrameExtractor.h>
#include <bob/ap/Energy.h>

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

#endif /* XBOB_AP_TYPES_H */
