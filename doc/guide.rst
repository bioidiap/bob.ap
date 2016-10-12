.. vim: set fileencoding=utf-8 :
.. Andre Anjos <andre.anjos@idiap.ch>
.. Pavel Korshunov <pavel.korshunov>@idiap.ch
.. Mon 17 Feb 2014 16:22:21 CET

.. testsetup:: aptest

  import os
  import sys
  import math
  import numpy
  import bob.ap

  import scipy.io.wavfile

  from pkg_resources import resource_filename
  wave_path = resource_filename('bob.ap', os.path.join('data', 'sample.wav'))

  sys.stdout =  open(os.devnull, 'w')
  rate, signal = scipy.io.wavfile.read(str(wave_path))
  sys.stdout = sys.__stdout__

************
 User Guide
************

This section will give more insight in simple and more complex
audio processing utilities of |project|. Currently, the following cepstral-based features are available:
using rectangular (RFCC), mel-scaled triangular (MFCC) [Davis1980]_, inverted mel-scaled triangular (IMFCC),
and linear triangular (LFCC) filters [Furui1981]_, spectral flux-based features (SSFC) [Scheirer1997]_,
subband centroid frequency (SCFC) [Le2011]_. We are planning to update and add more features in the
near future.

.. [Davis1980] S. Davis and P. Mermelstein, "Comparison of parametric representations for monosyllabic
   word recognition in continuously spoken sentences", in IEEE Transactions on Acoustics, Speech, and Signal Processing,
   1980, num 4, vol. 28, pages 357-366.
.. [Furui1981] S. Furui, Cepstral analysis technique for automatic speaker verification, in
   IEEE Transactions on Acoustics, Speech, and Signal Processing, 1981, num 2 vol 29, pages 254-272.
.. [Scheirer1997] E. Scheirer and M. Slaney, Construction and evaluation of a robust multifeature speech/music discriminator,
   in IEEE International Conference on Acoustics, Speech, and Signal Processing, ICASSP, 1997, vol 2, pages 1331-1334.
.. [Le2011] P. N. Le, E. Ambikairajah, J. Epps, V. Sethu, E. H. C. Choi, Investigation of Spectral Centroid Features for Cognitive Load Classification,
   in Speech Commun., April, 2011, num 4, vol 53, pages 540--551.

Simple audio processing
=======================

Below are 3 examples on how to read a wavefile and how to compute Linear frequency Cepstral Coefficients (LFCC) and Mel frequency cepstrum coefficients (MFCC).
Other features can be computed in a similar fashion (please check Python API for details).

Reading audio files
~~~~~~~~~~~~~~~~~~~

The usual native formats can be read with :py:mod:`scipy.io.wavfile` module.
These and other wave formats can be read through SoX_ using our native
bindings at :py:mod:`bob.io.audio`. An example of wave file can be found at
``bob/ap/test/data/sample.wav``.

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> import scipy.io.wavfile #doctest: +SKIP
  >>> rate, signal = scipy.io.wavfile.read(str(wave_path)) #doctest: +SKIP
  >>> print(rate)
  8000
  >>> print(signal)
  [  28   72   58 ..., -301   89  230]

In the above example, the sampling rate of the audio signal is **8 KHz** and
the signal array is of type **int16**.

User can directly compute the duration of signal (in seconds):

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> print(int(len(signal)/rate))
  2


LFCC and MFCC Extraction
~~~~~~~~~~~~~~~~~~~~~~~~

The LFCC and MFCC coefficients can be extracted from a audio signal by using
:py:class:`bob.ap.Ceps`. To do so, several parameters can be precised by the
user. Typically, these are precised in a configuration file. The following
values are the default ones:

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> win_length_ms = 20 # The window length of the cepstral analysis in milliseconds
  >>> win_shift_ms = 10 # The window shift of the cepstral analysis in milliseconds
  >>> n_filters = 24 # The number of filter bands
  >>> n_ceps = 19 # The number of cepstral coefficients
  >>> f_min = 0. # The minimal frequency of the filter bank
  >>> f_max = 4000. # The maximal frequency of the filter bank
  >>> delta_win = 2 # The integer delta value used for computing the first and second order derivatives
  >>> pre_emphasis_coef = 1.0 # The coefficient used for the pre-emphasis
  >>> dct_norm = True # A factor by which the cepstral coefficients are multiplied
  >>> mel_scale = True # Tell whether cepstral features are extracted on a linear (LFCC) or Mel (MFCC) scale

Once the parameters are precised, :py:class:`bob.ap.Ceps` can be called as
follows:

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> c = bob.ap.Ceps(rate, win_length_ms, win_shift_ms, n_filters, n_ceps, f_min, f_max, delta_win, pre_emphasis_coef, mel_scale, dct_norm)
  >>> signal = numpy.cast['float'](signal) # vector should be in **float**
  >>> mfcc = c(signal)
  >>> print(len(mfcc))
  199
  >>> print(len(mfcc[0]))
  19

LFCCs can be computed instead of MFCCs by setting ``mel_scale`` to ``False``:

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> c.mel_scale = False
  >>> lfcc = c(signal)

User can also choose to extract the energy. This is typically used for Voice
Activity Detection (VAD). Please check ``spkRecLib`` or ``FaceRecLib`` for more
details about VAD.

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> c.with_energy = True
  >>> lfcc_e = c(signal)
  >>> print(len(lfcc_e))
  199
  >>> print(len(lfcc_e[0]))
  20

It is also possible to compute first and second derivatives for those features:

.. doctest:: aptest
  :options: +NORMALIZE_WHITESPACE

  >>> c.with_delta = True
  >>> c.with_delta_delta = True
  >>> lfcc_e_d_dd = c(signal)
  >>> print(len(lfcc_e_d_dd))
  199
  >>> print(len(lfcc_e_d_dd[0]))
  60

.. include:: links.rst
