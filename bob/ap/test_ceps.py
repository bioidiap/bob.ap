#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Elie Khoury <Elie.Khoury@idiap.ch>
#
# Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland

import os
import numpy
import pkg_resources

from bob.sp import fft

from . import Ceps
from .test_utils import *

def cepstral_features_extraction(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta):

  #########################
  ## Initialisation part ##
  #########################
  c = Ceps(rate_wavsample[0], win_length_ms, win_shift_ms, n_filters, n_ceps, f_min, f_max, delta_win, pre_emphasis_coef)
  c.dct_norm = dct_norm
  c.mel_scale = mel_scale
  c.with_energy = with_energy
  c.with_delta = with_delta
  c.with_delta_delta = with_delta_delta
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
    #obj.assertAlmostEqual(ct.herz_to_mel(f_max), m_max, 7, "Error in Mel...")
    m_min = mel_python(f_min)
    #obj.assertAlmostEqual(ct.herz_to_mel(f_min), m_min, 7, "Error in Mel...")

    for i in range(n_filters + 2):
      alpha = ((i) / (n_filters + 1.0))
      f = mel_inv_python(m_min * (1 - alpha) + m_max * alpha)
      #obj.assertAlmostEqual(ct.mel_to_herz(m_min * (1 - alpha) + m_max * alpha), f, 7, "Error in MelInv...")
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
  ceps_sequence = numpy.zeros(n_ceps)
  dim0 = n_ceps
  if(with_energy):
    dim0 += + 1
  dim = dim0
  if(with_delta):
    dim += dim0
    if(with_delta_delta):
      dim += dim0
  else:
    with_delta_delta = False

  params = [ [ 0 for i in range(dim) ] for j in range(n_frames) ]

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

    if (with_energy):
      energy = sig_norm(win_length, frame, False)
      #e1 = ct.log_energy(frame)
      #obj.assertAlmostEqual(e1, energy, 7, "Error in Energy Computation...")

    f2 = numpy.copy(frame)

    # pre-emphasis filtering
    frame = pre_emphasis(frame, win_length, pre_emphasis_coef)
    #ct.pre_emphasis(f2)
    #for kk in range(len(frame)):
    #  obj.assertAlmostEqual(frame[kk], f2[kk], 7, "Error in Pre-Emphasis Computation...")

    # Hamming windowing
    f2 = numpy.copy(frame)
    frame = hamming_window(frame, hamming_kernel, win_length)
    #ct.hamming_window(f2)
    #for kk in range(len(frame)):
    #  obj.assertAlmostEqual(frame[kk], f2[kk], 7, "Error in Pre-Emphasis Computation...")

    f2=numpy.copy(frame)
    filters, x = log_filter_bank(frame, n_filters, p_index, win_size)

    #filt2 = ct.log_filter_bank(f2, win_size, n_filters)

    #for kk in range(len(filters)):
    #  obj.assertAlmostEqual(filters[kk], filt2[kk], 7, "Error in log Filtering")

    ceps = dct_transform(filters, n_filters, dct_kernel, n_ceps, dct_norm)
    #ceps2 = ct.apply_dct(n_ceps)


    if(with_energy):
      d1 = n_ceps + 1
      ceps[n_ceps] = energy
        #print(ceps)
    else:
      d1 = n_ceps

    # stock the results in params matrix
    vec=numpy.arange(d1)
    params[i][0:d1]=ceps[vec]

  # compute Delta coefficient
  if(with_delta):
    som = 0.0
    for i in range(1,delta_win+1):
      som = som + i*i
    som = som *2

    for i in range(n_frames):
      for k in range(n_ceps):
        params[i][d1+k] = 0.0
        for l in range(1, delta_win+1):
          if (i+l < n_frames):
            p_ind = i+l
          else:
            p_ind = n_frames - 1
          if (i-l > 0):
            n_ind = i-l
          else:
            n_ind = 0
          params[i][d1+k] = params[i][d1+k] + l * (params[p_ind][k] - params[n_ind][k])
        params[i][d1+k] = params[i][d1+k] / som

  # compute Delta of the Energy
  if(with_delta and with_energy):
    som = 0.0

    vec=numpy.arange(1,delta_win+1)
    som = 2.0* numpy.sum(vec*vec)

    for i in range(n_frames):
      k = n_ceps
      params[i][d1+k] = 0.0
      for l in range(1, delta_win+1):
        if (i+l < n_frames):
          p_ind = i+l
        else:
          p_ind = n_frames - 1
        if (i-l > 0):
          n_ind = i-l
        else:
          n_ind = 0
        params[i][d1+k] = params[i][d1+k] + l* (params[p_ind][k] - params[n_ind][k])
      params[i][d1+k] = params[i][d1+k] / som

  # compute Delta Delta of the coefficients
  if(with_delta_delta):
    som = 0.0
    for i in range(1,delta_win+1):
      som = som + i*i
    som = som *2
    for i in range(n_frames):
      for k in range(n_ceps):
        params[i][2*d1+k] = 0.0
        for l in range(1, delta_win+1):
          if (i+l < n_frames):
            p_ind = i+l
          else:
            p_ind = n_frames - 1
          if (i-l > 0):
            n_ind = i-l
          else:
            n_ind = 0
          params[i][2*d1+k] = params[i][2*d1+k] + l * (params[p_ind][d1+k] - params[n_ind][d1+k])
        params[i][2*d1+k] = params[i][2*d1+k] / som

  # compute Delta Delta of the energy
  if(with_delta_delta and with_energy):
    som = 0.0
    for i in range(1,delta_win+1):
      som = som + i*i
    som = som *2
    for i in range(n_frames):
      k = n_ceps
      params[i][2*d1+k] = 0.0
      for l in range(1, delta_win+1):
        if (i+l < n_frames):
          p_ind = i+l
        else:
          p_ind = n_frames - 1
        if (i-l > 0):
          n_ind = i-l
        else:
          n_ind = 0
        params[i][2*d1+k] = params[i][2*d1+k] + l * (params[p_ind][d1+k] - params[n_ind][d1+k])
      params[i][2*d1+k] = params[i][2*d1+k] / som
  data = numpy.array(params)

  return data

def cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                               pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta):
  c = Ceps(rate_wavsample[0], win_length_ms, win_shift_ms, n_filters, n_ceps, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, dct_norm)
  c.with_energy = with_energy
  c.with_delta = with_delta
  if c.with_delta:
    c.with_delta_delta = with_delta_delta
  #ct = TestCeps(c)
  A = c(rate_wavsample[1])
  B = cepstral_features_extraction(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm,
        f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)
  diff=numpy.sum(numpy.sum((A-B)*(A-B)))
  assert numpy.allclose(diff, 0., rtol=1e-07, atol=1e-05)

##################### Unit Tests ##################
def test_cepstral():
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

  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying win_length_ms
  win_length_ms = 30
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)
  # Varying win_shift_ms
  win_shift_ms = 15
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying n_filters
  n_filters = 20
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying n_ceps
  n_ceps = 12
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying f_min
  f_min = 300.
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying f_max
  f_max = 3300.
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying delta_win
  delta_win = 5
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying pre_emphasis_coef
  pre_emphasis_coef = 0.7
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying dct_norm
  dct_norm = False
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying mel_scale : linear scale
  mel_scale = False
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying with_energy
  with_energy = False
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying with_delta
  with_delta = False
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Varying with_delta_delta
  with_delta_delta = False
  cepstral_comparison_run(rate_wavsample, win_length_ms, win_shift_ms, n_filters, n_ceps, dct_norm, f_min, f_max, delta_win,
                             pre_emphasis_coef, mel_scale, with_energy, with_delta, with_delta_delta)

  # Test comparison operators and copy constructor
  c0 = Ceps(rate_wavsample[0], win_length_ms, win_shift_ms, n_filters, n_ceps, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, dct_norm)
  c1 = Ceps(c0)
  c2 = Ceps(c1)
  c2.win_length_ms = 27.
  assert c0 == c1
  assert not (c0 != c1)
  assert not (c0 == c2)
  assert c0 != c2
