#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Elie Khoury <Elie.Khoury@idiap.ch>
# Andre Anjos <andre.anjos@idiap.ch>
# Pavel Korshunov <Pavel.Korshunov@idiap.ch>
#
# Copyright (C) 2011-2014 Idiap Research Institute, Martigny, Switzerland

import os
import math
import numpy
import scipy.io.wavfile
import sys

from bob.sp import fft

def init_hamming_kernel(win_length):
  # Hamming initialisation
  cst = 2 * math.pi / (win_length - 1.0)
  hamming_kernel = numpy.zeros(win_length)

  for i in range(win_length):
    hamming_kernel[i] = (0.54 - 0.46 * math.cos(i * cst))
  return hamming_kernel

def init_freqfilter(rate, win_size,  mel_scale, n_filters, f_min, f_max):
  # Compute cut-off frequencies
  p_index = numpy.array(numpy.zeros(n_filters + 2), dtype=numpy.float64)
  if (mel_scale):
    # Mel scale
    m_max = mel_python(f_max)
    m_min = mel_python(f_min)

    for i in range(n_filters + 2):
      alpha = float(i) / (n_filters+1)
      f = mel_inv_python(m_min * (1 - alpha) + m_max * alpha)
      factor = float(f) / rate
      p_index[i] = win_size * factor
  else:
    # linear scale
    for i in range(n_filters + 2):
      alpha = float(i) / (n_filters+1)
      f = f_min * (1.0 - alpha) + f_max * alpha
      p_index[i] = float(win_size) / rate * f
  return p_index

def init_dct_kernel(n_filters, n_ceps, dct_norm):
  dct_kernel = numpy.zeros([n_ceps, n_filters], dtype=numpy.float64)

  dct_coeff = 1.0
  if dct_norm:
    dct_coeff = math.sqrt(2.0/n_filters)

  for i in range(0, n_ceps):
    for j in range(0, n_filters ):
      dct_kernel[i][j] = dct_coeff * math.cos(math.pi * i * (j + 0.5) / float(n_filters))

  if dct_norm:
    column_multiplier = numpy.ones(n_ceps, dtype=numpy.float64)
    column_multiplier[0] = math.sqrt(0.5)  # first element sqrt(0.5), the rest are 1.
    for j in range(0, n_filters):
      dct_kernel[:, j] = column_multiplier * dct_kernel[:, j]

  return dct_kernel

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

def pre_emphasis(frame, win_shift, coef, last_frame_elem):

  if (coef <= 0.0) or (coef > 1.0):
    print("Error: The emphasis coeff. should be between 0 and 1")
    return None

  last_element = frame[win_shift - 1]
  return numpy.append(frame[0]-coef * last_frame_elem, frame[1:]-coef*frame[:-1]), last_element

def hamming_window(vector, hamming_kernel, win_length):
  for i in range(win_length):
    vector[i] = vector[i] * hamming_kernel[i]
  return vector

def log_filter_bank(frame, n_filters, p_index, win_size):
  x1 = numpy.array(frame, dtype=numpy.complex128)
  complex_ = fft(x1)
  abscomplex = numpy.absolute(complex_)
  frame[0:int(win_size / 2) + 1] = abscomplex[0:int(win_size / 2) + 1]

  filters = log_triangular_bank(frame, n_filters, p_index)
  return filters, frame

def log_triangular_bank(data, n_filters, p_index):
  res_ = numpy.zeros(n_filters, dtype=numpy.float64)

  denominator = 1.0 / (p_index[1:n_filters+2] - p_index[0:n_filters+1])

  for i in range(0, n_filters):
    li = int(math.floor(p_index[i] + 1))
    mi = int(math.floor(p_index[i+1]))
    ri = int(math.floor(p_index[i+2]))
    if i == 0 or li == ri:
      li -= 1

    vec_left = numpy.arange(li, mi+1)
    vec_right = numpy.arange(mi+1, ri+1)
    res_[i] = numpy.sum(data[vec_left] * denominator[i] * (vec_left-p_index[i])) + \
              numpy.sum(data[vec_right] * denominator[i+1] * (p_index[i+2]-vec_right))
    # alternative but equivalent implementation:
    # filt = numpy.zeros(ri-li+1, dtype=numpy.float64)
    # filt_l = denominator[i] * (vec_left-p_index[i])
    # filt_p = denominator[i+1] * (p_index[i+2]-vec_right)
    # filt = numpy.append(filt_l, filt_p)
    # vect_full = numpy.arange(li, ri+1)
    # res_[i] = numpy.sum(data[vect_full] * filt)

  FBANK_OUT_FLOOR = sys.float_info.epsilon
  return numpy.log(numpy.where(res_ < FBANK_OUT_FLOOR, FBANK_OUT_FLOOR, res_))

def dct_transform(filters, n_filters, dct_kernel, n_ceps):

  ceps = numpy.zeros(n_ceps)
  vec = numpy.array(range(0, n_filters))
  for i in range(0, n_ceps):
    ceps[i] = numpy.sum(filters[vec] * dct_kernel[i])

  return ceps
