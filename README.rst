.. vim: set fileencoding=utf-8 :
.. Andre Anjos <andre.anjos@idiap.ch>
.. Thu 30 Jan 08:46:53 2014 CET

=============================
 Python bindings for bob.ap
=============================

This package contains a set of Pythonic bindings for Bob's signal processing
package and functionality.

Installation
------------

Install it through normal means, via PyPI or use ``zc.buildout`` to bootstrap
the package and run test units.

Documentation
-------------

You can generate the documentation for this package, after installation, using
Sphinx::

  $ aphinx-build -b html doc aphinx

This shall place in the directory ``aphinx``, the current version for the
documentation of the package.

Testing
-------

You can run a set of tests using the nose test runner::

  $ nosetests -sv xbob.ap

.. warning::

   If Bob <= 1.2.1 is installed on your python path, nose will automatically
   load the old version of the insulate plugin available in Bob, which will
   trigger the loading of incompatible shared libraries (from Bob itself), in
   to your working binary. This will cause a stack corruption. Either remove
   the centrally installed version of Bob, or build your own version of Python
   in which Bob <= 1.2.1 is not installed.

You can run our documentation tests using aphinx itself::

  $ aphinx-build -b doctest doc aphinx

You can test overall test coverage with::

  $ nosetests --with-coverage --cover-package=xbob.ap

The ``coverage`` egg must be installed for this to work properly.

Development
-----------

To develop this package, install using ``zc.buildout``, using the buildout
configuration found on the root of the package::

  $ python bootstrap.py
  ...
  $ ./bin/buildout

Tweak the options in ``buildout.cfg`` to disable/enable verbosity and debug
builds.
