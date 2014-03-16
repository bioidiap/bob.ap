#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Elie Khoury <Elie.Khoury@idiap.ch>
#
# Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland

import os
import numpy
import math
import pkg_resources

from . import Energy
from .test_utils import *

def energy_computation(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta):

  #########################
  ## Initialisation part ##
  #########################

  c = Energy(rate_wavsample[0], win_length_ms, win_shift_ms)

  #ct = TestCeps(c)

  sf = rate_wavsample[0]
  data = rate_wavsample[1]

  win_length = int (sf * win_length_ms / 1000)
  win_shift = int (sf * win_shift_ms / 1000)
  win_size = int (2.0 ** math.ceil(math.log(win_length) / math.log(2)))
  m = int (math.log(win_size) / math.log(2))

  # Hamming initialisation
  cst = 2 * math.pi / (win_length - 1.0)
  hamming_kernel = numpy.zeros(win_length)

  for i in range(win_length):
    hamming_kernel[i] = (0.54 - 0.46 * math.cos(i * cst))

  # Compute cut-off frequencies
  p_index = numpy.array(numpy.zeros(n_filters + 2), dtype=numpy.int16)
  if(mel_scale):
    # Mel scale
    m_max = mel_python(f_max)
    m_min = mel_python(f_min)

    for i in range(n_filters + 2):
      alpha = ((i) / (n_filters + 1.0))
      f = mel_inv_python(m_min * (1 - alpha) + m_max * alpha)
      factor = f / (sf * 1.0)
      p_index[i] = int (round((win_size) * factor))
  else:
    #linear scale
    for i in range(n_filters + 2):
      alpha = (i) / (n_filters + 1.0)
      f = f_min * (1.0 - alpha) + f_max * alpha
      p_index[i] = int (round((win_size / (sf * 1.0) * f)))

  #Cosine transform initialisation
  dct_kernel = [ [ 0 for i in range(n_filters) ] for j in range(n_ceps) ]

  for i in range(1, n_ceps + 1):
    for j in range(1, n_filters + 1):
      dct_kernel[i - 1][j - 1] = math.cos(math.pi * i * (j - 0.5) / n_filters)

  ######################################
  ### End of the Initialisation part ###
  ######################################

  ######################################
  ###          Core code             ###
  ######################################

  data_size = data.shape[0]
  n_frames = int(1 + (data_size - win_length) / win_shift)

  # create features set

  params = [ 0 for j in range(n_frames) ]

  # compute cepstral coefficients
  delta = 0
  for i in range(n_frames):
    # create a frame
    frame = numpy.zeros(win_size, dtype=numpy.float64)
    som = 0.0
    vec = numpy.arange(win_length)
    frame[vec] = data[vec + i * win_shift]
    som = numpy.sum(frame)
    som = som / win_size
    frame = frame - som


    energy = sig_norm(win_length, frame, False)

    params[i] = energy

  data = numpy.array(params)

  return data


def energy_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta):
  c = Energy(rate_wavsample[0], win_length_ms, win_shift_ms)

  #ct = TestCeps(c)
  A = c(rate_wavsample[1])

  B = energy_computation(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm,
        f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  diff=numpy.sum(numpy.sum((A-B)*(A-B)))
  assert numpy.allclose(diff, 0., rtol=1e-07, atol=1e-05)

##################### Unit Tests ##################
def test_energy():

  rate_wavsample = read(pkg_resources.resource_filename(__name__, os.path.join('data', 'sample.wav')))

  win_length_ms = 20
  win_shift_ms = 10
  n_filters = 24
  n_ceps = 19
  f_min = 0.
  f_max = 4000.
  delta_win = 2
  pre_emphasis_coef = 0.97
  dct_norm = True
  mel_scale = True
  with_energy = True
  with_delta = True
  with_delta_delta = True

  energy_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)
