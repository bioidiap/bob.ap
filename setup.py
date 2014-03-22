#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Andre Anjos <andre.anjos@idiap.ch>
# Thu 30 Jan 08:45:49 2014 CET

from setuptools import setup, find_packages, dist
dist.Distribution(dict(setup_requires=['xbob.blitz']))
from xbob.blitz.extension import Extension

packages = ['bob-ap >= 1.2.2']
version = '2.0.0a0'

setup(

    name='xbob.ap',
    version=version,
    description='Bindings for Bob\'s audio processing utilities',
    url='http://github.com/anjos/xbob.ap',
    license='BSD',
    author='Andre Anjos',
    author_email='andre.anjos@idiap.ch',

    long_description=open('README.rst').read(),

    packages=find_packages(),
    include_package_data=True,

    install_requires=[
      'setuptools',
      'xbob.blitz',
      'xbob.sp', # for testing
      'scipy', # for testing
    ],

    namespace_packages=[
      "xbob",
      ],

    ext_modules = [
      Extension("xbob.ap.version",
        [
          "xbob/ap/version.cpp",
          ],
        packages = packages,
        version = version,
        ),
      Extension("xbob.ap._library",
        [
          "xbob/ap/energy.cpp",
          "xbob/ap/frame_extractor.cpp",
          "xbob/ap/spectrogram.cpp",
          "xbob/ap/ceps.cpp",
          "xbob/ap/main.cpp",
          ],
        packages = packages,
        version = version,
        ),
      ],

    classifiers = [
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: BSD License',
      'Natural Language :: English',
      'Programming Language :: Python',
      'Programming Language :: Python :: 3',
      'Topic :: Software Development :: Libraries :: Python Modules',
      ],

    )
