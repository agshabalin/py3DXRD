# -*- coding: utf-8 -*-
import os, sys, setuptools
from distutils.core import setup
import numpy as np

setup(name='Geometry Converter',
      version=get_version(),
      author='Anatoly Shabalin',
      author_email='a.shabalin.r@egmail.com',
      cmdclass={ 'build_ext' : build_ext_subclass },
      description='Converts geometry between HEXRD, ImageD11, PolyXSim, and pyFAI.',
      long_description = open('README.md').read(),
      license = "GPL",
      ext_package = "ImageD11",
      setup_requires = ["numpy", "scipy", "setuptools"],   # to compile
      install_requires = ["numpy", "scipy", "setuptools"],
      extras_require = { 'full' : [], 'rare' : [] },
      packages = [],
      package_dir = {"geometry_converter":"geometry_converter"},
      url = "https://github.com/agshabalin/py3DXRD.git",
      scripts = ["py3DXRD/convert_geometry.py"],
      long_description_content_type='text/markdown',
)