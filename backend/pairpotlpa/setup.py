#Building LPA.pyd file
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import pkg_resources

ext_modules = [

    Extension('pairpotlpa',
              sources=['LPA.cpp'],
              language='c++',
              extra_compile_args=['-std=c++11'],
              include_dirs = [pkg_resources.resource_filename('pybind11', 'include/')],
              extra_link_args=[]),
]

setup(name='pairpotlpa',
      version='0.1',
      author="rzh, zzj",
      author_email="rrrzhan@mail.nankai.edu.cn",
      description=("Cell Propagation algorithm implemented in C++"),
      ext_modules=ext_modules,
      cmdclass={'build_ext': build_ext},
      zip_safe=False)
