#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Elie Khoury <Elie.Khoury@idiap.ch>
# Andre Anjos <andre.anjos@idiap.ch>
#
# Copyright (C) 2011-2014 Idiap Research Institute, Martigny, Switzerland

import os
import math
import numpy
import scipy.io.wavfile

from bob.sp import fft

def read(filename):
  """Read video.FrameContainer containing preprocessed frames"""

  fileName, fileExtension = os.path.splitext(filename)
  wav_filename = filename
  rate, data = scipy.io.wavfile.read(str(wav_filename)) # the data is read in its native format
  if data.dtype =='int16':
    data = numpy.cast['float'](data)
  return [rate,data]

def compare(v1, v2, width):
  return abs(v1-v2) <= width

def mel_python(f):
  return 2595.0*math.log10(1.+f/700.0)

def mel_inv_python(value):
  return 700.0 * (10 ** (value / 2595.0) - 1)

def sig_norm(win_length, frame, flag):
  gain = 0.0
  for i in range(win_length):
    gain = gain + frame[i] * frame[i]

  ENERGY_FLOOR = 1.0
  if gain < ENERGY_FLOOR:
    gain = math.log(ENERGY_FLOOR)
  else:
    gain = math.log(gain)

  if(flag and gain != 0.0):
    for i in range(win_length):
      frame[i] = frame[i] / gain
  return gain

def pre_emphasis(frame, win_length, a):
  if (a < 0.0) or (a >= 1.0):
    print("Error: The emphasis coeff. should be between 0 and 1")
  if (a == 0.0):
    return frame
  else:
    for i in range(win_length - 1, 0, -1):
      frame[i] = frame[i] - a * frame[i - 1]
    frame[0] = (1. - a) * frame[0]
  return frame

def hamming_window(vector, hamming_kernel, win_length):
  for i in range(win_length):
    vector[i] = vector[i] * hamming_kernel[i]
  return vector

def log_filter_bank(x, n_filters, p_index, win_size):

  x1 = numpy.array(x, dtype=numpy.complex128)
  complex_ = fft(x1)
  for i in range(0, int(win_size / 2) + 1):
    re = complex_[i].real
    im = complex_[i].imag
    x[i] = math.sqrt(re * re + im * im)
  filters = log_triangular_bank(x, n_filters, p_index)
  return filters, x

def log_triangular_bank(data, n_filters, p_index):
  a = 1.0 / (p_index[1:n_filters+2] - p_index[0:n_filters+1] + 1)
  vec1 =  list(numpy.arange(p_index[i], p_index[i + 1]) for i in range(0, n_filters))
  vec2 =  list(numpy.arange(p_index[i+1], p_index[i + 2] + 1) for i in range(0, n_filters))
  res_ = numpy.array([(numpy.sum(data[vec1[i]]*(1.0 - a [i]* (p_index[i + 1]-(vec1[i])))) +
          numpy.sum(data[vec2[i]] * (1.0 - a[i+1] * ( (vec2[i]) - p_index[i + 1]))))
          for i in range(0, n_filters)])
  FBANK_OUT_FLOOR = 1.0
  filters = numpy.log(numpy.where(res_ < FBANK_OUT_FLOOR, FBANK_OUT_FLOOR, res_))
  return filters

def dct_transform(filters, n_filters, dct_kernel, n_ceps, dct_norm):
  if dct_norm:
    dct_coeff = numpy.sqrt(2.0/(n_filters))
  else :
    dct_coeff = 1.0

  ceps = numpy.zeros(n_ceps + 1)
  vec = numpy.array(range(1, n_filters + 1))
  for i in range(1, n_ceps + 1):
    ceps[i - 1] = numpy.sum(filters[vec - 1] * dct_kernel[i - 1][0:n_filters])
    ceps[i - 1] = ceps[i - 1] * dct_coeff

  return ceps
