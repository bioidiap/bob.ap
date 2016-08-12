#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
# Andre Anjos <andre.anjos@idiap.ch>
# Thu 30 Jan 08:45:49 2014 CET

bob_packages = ['bob.core', 'bob.sp']

from setuptools import setup, find_packages, dist
dist.Distribution(dict(setup_requires=['bob.extension', 'bob.blitz'] + bob_packages))
from bob.blitz.extension import Extension, Library, build_ext

from bob.extension.utils import load_requirements
build_requires = load_requirements()

version = open('version.txt').read().rstrip()

setup(

    name='bob.ap',
    version=version,
    description='Bob\'s audio processing utilities',
    url='https://gitlab.idiap.ch/bob/bob.ap',
    license='BSD',
    author='Andre Anjos',
    author_email='andre.anjos@idiap.ch',

    long_description=open('README.rst').read(),

    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,

    setup_requires = build_requires,
    install_requires = build_requires,

    ext_modules = [
      Extension("bob.ap.version",
        [
          "bob/ap/version.cpp",
        ],
        version = version,
        bob_packages = bob_packages,
      ),

      Library("bob.ap.bob_ap",
        [
          "bob/ap/cpp/Energy.cpp",
          "bob/ap/cpp/FrameExtractor.cpp",
          "bob/ap/cpp/Spectrogram.cpp",
          "bob/ap/cpp/Ceps.cpp",
        ],
        version = version,
        bob_packages = bob_packages,
      ),

      Extension("bob.ap._library",
        [
          "bob/ap/energy.cpp",
          "bob/ap/frame_extractor.cpp",
          "bob/ap/spectrogram.cpp",
          "bob/ap/ceps.cpp",
          "bob/ap/main.cpp",
          ],
        version = version,
        bob_packages = bob_packages,
      ),
    ],

    cmdclass = {
      'build_ext': build_ext
    },

    classifiers = [
      'Framework :: Bob',
      'Development Status :: 4 - Beta',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: BSD License',
      'Natural Language :: English',
      'Programming Language :: Python',
      'Programming Language :: Python :: 3',
      'Topic :: Software Development :: Libraries :: Python Modules',
    ],

  )
