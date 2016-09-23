.. vim: set fileencoding=utf-8 :
.. Andre Anjos <andre.anjos@idiap.ch>
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

This section will give a deeper insight in some simple and some more complex
audio processing utilities of |project|. Currently, only cepstral extraction
module is available. We are planning to update and add more features in the
near future.

Simple audio processing
=======================

Below are 3 examples on how to read a wavefile and how to compute Linear frequency Cepstral Coefficients (LFCC) and Mel frequency cepstrum coefficients (MFCC).

Reading audio files
~~~~~~~~~~~~~~~~~~~~

The usual native formats can be read with :py:mod:`scipy.io.wavfile` module. Other
wave formats can be found in some other python modules like :py:mod:`pysox`. An
example of wave file can be found here ``bob/ap/test/data/sample.wav``

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
  >>> n_filters = 20 # The number of filter bands
  >>> n_ceps = 20 # The number of cepstral coefficients
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

