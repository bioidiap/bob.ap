.. vim: set fileencoding=utf-8 :
.. Sat 13 Aug 2016 00:24:52 CEST

.. image:: http://img.shields.io/badge/docs-stable-yellow.png
   :target: http://pythonhosted.org/bob.ap/index.html
.. image:: http://img.shields.io/badge/docs-latest-orange.png
   :target: https://www.idiap.ch/software/bob/docs/latest/bob/bob.ap/master/index.html
.. image:: https://gitlab.idiap.ch/bob/bob.ap/badges/v2.1.1/build.svg
   :target: https://gitlab.idiap.ch/bob/bob.ap/commits/v2.1.1
.. image:: https://img.shields.io/badge/gitlab-project-0000c0.svg
   :target: https://gitlab.idiap.ch/bob/bob.ap
.. image:: http://img.shields.io/pypi/v/bob.ap.png
   :target: https://pypi.python.org/pypi/bob.ap
.. image:: http://img.shields.io/pypi/dm/bob.ap.png
   :target: https://pypi.python.org/pypi/bob.ap


==========================
 Audio Processing for Bob
==========================

This package is part of the signal-processing and machine learning toolbox
Bob_. It contains basic audio processing utilities. Currently, the following cepstral-based features are available:
using rectangular (RFCC), mel-scaled triangular (MFCC) [Davis1980]_, inverted mel-scaled triangular (IMFCC),
and linear triangular (LFCC) filters [Furui1981]_, spectral flux-based features (SSFC) [Scheirer1997]_,
subband centroid frequency (SCFC) [Le2011]_. We are planning to update and add more features in the
near future.

*Please note that the implementation of MFCC and LFCC features has changed compared to an earlier version of the package,
as we corrected pre-emphasis and DCT computations. Delta and delta-delta computations were slightly changed too.*

Installation
------------

Follow our `installation`_ instructions. Then, using the Python interpreter
provided by the distribution, bootstrap and buildout this package::

  $ python bootstrap-buildout.py
  $ ./bin/buildout


Contact
-------

For questions or reporting issues to this software package, contact our
development `mailing list`_.

.. [Davis1980] S. Davis and P. Mermelstein, "Comparison of parametric representations for monosyllabic
   word recognition in continuously spoken sentences", in IEEE Transactions on Acoustics, Speech, and Signal Processing,
   1980, num 4, vol. 28, pages 357-366.
.. [Furui1981] S. Furui, Cepstral analysis technique for automatic speaker verification, in
   IEEE Transactions on Acoustics, Speech, and Signal Processing, 1981, num 2 vol 29, pages 254-272.
.. [Scheirer1997] E. Scheirer and M. Slaney, Construction and evaluation of a robust multifeature speech/music discriminator,
   in IEEE International Conference on Acoustics, Speech, and Signal Processing, ICASSP, 1997, vol 2, pages 1331-1334.
.. [Le2011] P. N. Le, E. Ambikairajah, J. Epps, V. Sethu, E. H. C. Choi, Investigation of Spectral Centroid Features for Cognitive Load Classification,
   in Speech Commun., April, 2011, num 4, vol 53, pages 540--551.

.. Place your references here:
.. _bob: https://www.idiap.ch/software/bob
.. _installation: https://gitlab.idiap.ch/bob/bob/wikis/Installation
.. _mailing list: https://groups.google.com/forum/?fromgroups#!forum/bob-devel

