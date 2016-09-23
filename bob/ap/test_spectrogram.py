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

from scipy import signal

from . import Spectrogram
from .test_utils import *

def spectrogram_computation(rate_wavsample, win_length_ms, win_shift_ms, n_filters,
                            f_min, f_max, pre_emphasis_coef, mel_scale, normalize_mean):
  #########################
  ## Initialisation part ##
  #########################

  rate = rate_wavsample[0]
  data = rate_wavsample[1]

  win_length = int (rate * win_length_ms / 1000)
  win_shift = int (rate * win_shift_ms / 1000)
  win_size = int (2.0 ** math.ceil(math.log(win_length) / math.log(2)))
  m = int (math.log(win_size) / math.log(2))

  # Hamming initialisation
  hamming_kernel = init_hamming_kernel(win_length)

  # Compute cut-off frequencies
  p_index = init_freqfilter(rate, win_size,  mel_scale, n_filters, f_min, f_max)

  ######################################
  ### End of the Initialisation part ###
  ######################################

  ######################################
  ###          Core code             ###
  ######################################

  data_size = data.shape[0]
  n_frames = int(1 + (data_size - win_length) / win_shift)

  # create features set
  features = numpy.zeros([n_frames, int(win_size/2)+1], dtype=numpy.float64)

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

    frame_, last_frame_elem = pre_emphasis(frame[vec], win_shift, pre_emphasis_coef, last_frame_elem)
    frame[vec] = frame_

    # Hamming windowing
    frame = hamming_window(frame, hamming_kernel, win_length)

    filters, spec_row = log_filter_bank(frame, n_filters, p_index, win_size)

    features[i] = spec_row[0:int(win_size/2)+1]

  return numpy.array(features)

def spectrogram_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, f_min, f_max,
                               pre_emphasis_coef, mel_scale, normalize_mean):
  c = Spectrogram(rate_wavsample[0], win_length_ms, win_shift_ms, n_filters,
                  f_min, f_max, pre_emphasis_coef, mel_scale, normalize_mean)

  A = c(rate_wavsample[1])
  B = spectrogram_computation(rate_wavsample, win_length_ms, win_shift_ms, n_filters,
                              f_min, f_max, pre_emphasis_coef, mel_scale, normalize_mean)

  diff=numpy.sum(numpy.sum((A-B)*(A-B)))
  assert numpy.allclose(diff, 0., rtol=1e-07, atol=1e-05)


##################### Unit Tests ##################
def test_spectrogram():

    rate_wavsample = read(pkg_resources.resource_filename(__name__, os.path.join('data', 'sample.wav')))

    win_length_ms = 20
    win_shift_ms = 10
    normalize_mean = True
    n_filters = 20
    f_min = 0.
    f_max = 4000.
    pre_emphasis_coef = 1.0
    mel_scale = True
    spectrogram_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, f_min, f_max,
                               pre_emphasis_coef, mel_scale, normalize_mean)
