#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Elie Khoury <Elie.Khoury@idiap.ch>
# Pavel Korshunov <Pavel.Korshunov@idiap.ch>
#
# Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland

import os
import numpy
import pkg_resources

from bob.sp import fft

from . import Ceps
from .test_utils import *

class DotDict(dict):
    def __getattr__(self, key):
        return self[key]
    def __setattr__(self, key, val):
        if key in self.__dict__:
            self.__dict__[key] = val
        else:
            self[key] = val


cepstral_params = DotDict(dict(win_length_ms=20,
                               win_shift_ms=10,
                               n_filters=20,
                               f_min=0.,
                               f_max=4000.,
                               pre_emphasis_coef=1.0,
                               mel_scale=True,
                               n_ceps=20,
                               delta_win=2,
                               dct_norm=True,
                               normalize_mean=True,
                               with_energy=True,
                               with_delta=True,
                               with_delta_delta=True
                               ))


def cepstral_python_implementation(rate_wavsample, cepstral_params):

  #########################
  ## Initialisation part ##
  #########################
  rate = rate_wavsample[0]
  data = rate_wavsample[1]

  win_length = int(rate * cepstral_params.win_length_ms / 1000)
  win_shift = int(rate * cepstral_params.win_shift_ms / 1000)
  win_size = int(2.0 ** math.ceil(math.log(win_length) / math.log(2)))
  m = int(math.log(win_size) / math.log(2))

  # Hamming initialisation
  hamming_kernel = init_hamming_kernel(win_length)

  # Compute cut-off frequencies
  p_index = init_freqfilter(rate, win_size, cepstral_params.mel_scale, cepstral_params.n_filters,
                            cepstral_params.f_min, cepstral_params.f_max)

  #Cosine transform initialisation
  dct_kernel = init_dct_kernel(cepstral_params.n_filters, cepstral_params.n_ceps, cepstral_params.dct_norm)

  ######################################
  ### End of the Initialisation part ###
  ######################################

  ######################################
  ###          Core code             ###
  ######################################

  data_size = data.shape[0]
  n_frames = int(1 + (data_size - win_length) / win_shift)

  # create features set
  dim0 = cepstral_params.n_ceps
  if(cepstral_params.with_energy):
    dim0 += + 1
  dim = dim0
  if(cepstral_params.with_delta):
    dim += dim0
    if(cepstral_params.with_delta_delta):
      dim += dim0
  else:
    cepstral_params.with_delta_delta = False

  features = numpy.zeros([n_frames, dim], dtype=numpy.float64)

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

    if (cepstral_params.with_energy):
      energy = sig_norm(win_length, frame, False)

    # pre-emphasis filtering
    frame_, last_frame_elem = pre_emphasis(frame[vec], win_shift, cepstral_params.pre_emphasis_coef, last_frame_elem)
    frame[vec] = frame_

    # Hamming windowing
    frame = hamming_window(frame, hamming_kernel, win_length)

    # FFT and filters
    filters, spec_row = log_filter_bank(frame, cepstral_params.n_filters, p_index, win_size)

    # apply DCT
    ceps = dct_transform(filters, cepstral_params.n_filters, dct_kernel, cepstral_params.n_ceps)


    ######################################
    ###     Deltas and Delta-Deltas    ###
    ######################################

    d1 = cepstral_params.n_ceps
    if cepstral_params.with_energy:
      d1 = cepstral_params.n_ceps + 1
      ceps = numpy.append(ceps, energy)

    # stock the results in features matrix
    vec=numpy.arange(d1)
    features[i][0:d1]=ceps[vec]

  # compute Delta coefficient
  if(cepstral_params.with_delta):
    som = 0.0
    for i in range(1,cepstral_params.delta_win+1):
      som = som + i*i
    som = som *2

    for i in range(n_frames):
      for k in range(cepstral_params.n_ceps):
        features[i][d1+k] = 0.0
        for l in range(1, cepstral_params.delta_win+1):
          if (i+l < n_frames):
            p_ind = i+l
          else:
            p_ind = n_frames - 1
          if (i-l > 0):
            n_ind = i-l
          else:
            n_ind = 0
          features[i][d1+k] = features[i][d1+k] + l * (features[p_ind][k] - features[n_ind][k])
        # features[i][d1+k] = features[i][d1+k] / som  # do not normalize anymore

  # compute Delta of the Energy
  if(cepstral_params.with_delta and cepstral_params.with_energy):
    som = 0.0

    vec=numpy.arange(1,cepstral_params.delta_win+1)
    som = 2.0* numpy.sum(vec*vec)

    for i in range(n_frames):
      k = cepstral_params.n_ceps
      features[i][d1+k] = 0.0
      for l in range(1, cepstral_params.delta_win+1):
        if (i+l < n_frames):
          p_ind = i+l
        else:
          p_ind = n_frames - 1
        if (i-l > 0):
          n_ind = i-l
        else:
          n_ind = 0
        features[i][d1+k] = features[i][d1+k] + l* (features[p_ind][k] - features[n_ind][k])
      # features[i][d1+k] = features[i][d1+k] / som  # do not normalize anymore

  # compute Delta Delta of the coefficients
  if(cepstral_params.with_delta_delta):
    som = 0.0
    for i in range(1,cepstral_params.delta_win+1):
      som = som + i*i
    som = som *2
    for i in range(n_frames):
      for k in range(cepstral_params.n_ceps):
        features[i][2*d1+k] = 0.0
        for l in range(1, cepstral_params.delta_win+1):
          if (i+l < n_frames):
            p_ind = i+l
          else:
            p_ind = n_frames - 1
          if (i-l > 0):
            n_ind = i-l
          else:
            n_ind = 0
          features[i][2*d1+k] = features[i][2*d1+k] + l * (features[p_ind][d1+k] - features[n_ind][d1+k])
        # features[i][2*d1+k] = features[i][2*d1+k] / som  # do not normalize anymore

  # compute Delta Delta of the energy
  if(cepstral_params.with_delta_delta and cepstral_params.with_energy):
    som = 0.0
    for i in range(1,cepstral_params.delta_win+1):
      som = som + i*i
    som = som *2
    for i in range(n_frames):
      k = cepstral_params.n_ceps
      features[i][2*d1+k] = 0.0
      for l in range(1, cepstral_params.delta_win+1):
        if (i+l < n_frames):
          p_ind = i+l
        else:
          p_ind = n_frames - 1
        if (i-l > 0):
          n_ind = i-l
        else:
          n_ind = 0
        features[i][2*d1+k] = features[i][2*d1+k] + l * (features[p_ind][d1+k] - features[n_ind][d1+k])
      # features[i][2*d1+k] = features[i][2*d1+k] / som  # do not normalize anymore

  return numpy.array(features)

def cepstral_comparison_run(cepstral_params):
  rate_wavsample = read(pkg_resources.resource_filename(__name__, os.path.join('data', 'sample.wav')))

  cepstral = Ceps(rate_wavsample[0], cepstral_params.win_length_ms, cepstral_params.win_shift_ms,
                  cepstral_params.n_filters, cepstral_params.n_ceps, cepstral_params.f_min, cepstral_params.f_max,
                  cepstral_params.delta_win, cepstral_params.pre_emphasis_coef, cepstral_params.mel_scale,
                  cepstral_params.dct_norm, cepstral_params.normalize_mean)
  cepstral.with_energy = cepstral_params.with_energy
  cepstral.with_delta = cepstral_params.with_delta
  if cepstral.with_delta:
    cepstral.with_delta_delta = cepstral_params.with_delta_delta

  A = cepstral(rate_wavsample[1])
  B = cepstral_python_implementation(rate_wavsample, cepstral_params)

  diff=numpy.sum(numpy.sum((A-B)*(A-B)))
  assert numpy.allclose(diff, 0., rtol=1e-07, atol=1e-05)

##################### Unit Tests ##################
def test_mfcc1():
  cepstral_comparison_run(cepstral_params)

def test_mfcc2():
  # Varying win_length_ms
  cepstral_params.win_length_ms = 30
  cepstral_comparison_run(cepstral_params)

def test_mfcc3():
  # Varying win_shift_ms
  cepstral_params.win_shift_ms = 15
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc4():
  # Varying n_filters
  cepstral_params.n_filters = 24
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc5():
  # Varying n_ceps
  cepstral_params.n_ceps = 19
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc6():
  # Varying f_min
  cepstral_params.f_min = 300.
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc7():
  # Varying f_max
  cepstral_params.f_max = 3300.
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc8():
  # Varying delta_win
  cepstral_params.delta_win = 5
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc9():
  # Varying pre_emphasis_coef
  cepstral_params.pre_emphasis_coef = 0.95
  cepstral_comparison_run(cepstral_params)
                             

def test_mfcc10():
  # Varying dct_norm
  cepstral_params.dct_norm = False
  cepstral_comparison_run(cepstral_params)
                             

def test_lfcc1():
  # Varying mel_scale : linear scale
  cepstral_params.mel_scale = False
  cepstral_comparison_run(cepstral_params)
                             

def test_lfcc2():
  # Varying with_energy
  cepstral_params.with_energy = False
  cepstral_comparison_run(cepstral_params)
                             

def test_lfcc3():
  # Varying with_delta
  cepstral_params.with_delta_delta = False
  cepstral_comparison_run(cepstral_params)
                             

def test_lfcc4():
  # Varying with_delta_delta
  cepstral_params.with_delta = False
  cepstral_comparison_run(cepstral_params)

def test_cepstral_copy():
  # # Test comparison operators and copy constructor
  rate_wavsample = read(pkg_resources.resource_filename(__name__, os.path.join('data', 'sample.wav')))
  cepstral = Ceps(rate_wavsample[0], cepstral_params.win_length_ms, cepstral_params.win_shift_ms,
                  cepstral_params.n_filters, cepstral_params.n_ceps, cepstral_params.f_min,
                  cepstral_params.f_max, cepstral_params.delta_win, cepstral_params.pre_emphasis_coef,
                  cepstral_params.mel_scale, cepstral_params.dct_norm, cepstral_params.normalize_mean, )
  c1 = Ceps(cepstral)
  c2 = Ceps(c1)
  c2.win_length_ms = 27.
  assert cepstral == c1
  assert not (cepstral != c1)
  assert not (cepstral == c2)
  assert cepstral != c2
