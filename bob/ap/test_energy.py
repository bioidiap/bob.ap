#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Elie Khoury <Elie.Khoury@idiap.ch>
# Pavel Korshunov <Pavel.Korshunov@idiap.ch>
#
# Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland

import os
import numpy
import math
import pkg_resources

from . import Energy
from .test_utils import *

def energy_computation(rate_wavsample, win_length_ms, win_shift_ms, normalize_mean):

  #########################
  ## Initialisation part ##
  #########################

  rate = rate_wavsample[0]
  data = rate_wavsample[1]

  win_length = int (rate * win_length_ms / 1000)
  win_shift = int (rate * win_shift_ms / 1000)
  win_size = int (2.0 ** math.ceil(math.log(win_length) / math.log(2)))

  ######################################
  ### End of the Initialisation part ###
  ######################################

  ######################################
  ###          Core code             ###
  ######################################

  data_size = data.shape[0]
  n_frames = int(1 + (data_size - win_length) / win_shift)

  # create features set

  features = [ 0 for j in range(n_frames) ]

  last_frame_elem = 0
  # compute cepstral coefficients
  for i in range(n_frames):
    # create a frame
    frame = numpy.zeros(win_size, dtype=numpy.float64)
    vec = numpy.arange(win_length)
    frame[vec] = data[vec + i * win_shift]
    som = numpy.sum(frame)
    som = som / win_size
    frame[vec] -= som  # normalization by mean here

    energy = sig_norm(win_length, frame, False)
    features[i] = energy

  return numpy.array(features)


def energy_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, normalize_mean):
  c = Energy(rate_wavsample[0], win_length_ms, win_shift_ms, normalize_mean)

  A = c(rate_wavsample[1])

  B = energy_computation(rate_wavsample, win_length_ms, win_shift_ms, normalize_mean)

  diff=numpy.sum(numpy.sum((A-B)*(A-B)))
  assert numpy.allclose(diff, 0., rtol=1e-07, atol=1e-05)

##################### Unit Tests ##################
def test_energy():

  rate_wavsample = read(pkg_resources.resource_filename(__name__, os.path.join('data', 'sample.wav')))

  win_length_ms = 20
  win_shift_ms = 10
  normalize_mean = True

  energy_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, normalize_mean)
