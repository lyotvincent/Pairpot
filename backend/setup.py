#Building LPA.pyd file
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os

ext_modules = [

    Extension('label_propagation',
              sources=['LPA.cpp'],
              language='c++',
              extra_compile_args=['-std=c++11'],
              include_dirs=['/home/rzh/Browser/lib/python3.10/site-packages/pybind11/include/'],
              extra_link_args=[]),
]

setup(name='label_propagation',
      version='0.1',
      author="rzh, zzj",
      author_email="rrrzhan@mail.nankai.edu.cn",
      description=("A label propagation algorithm implemented in C++"),
      ext_modules=ext_modules,
      cmdclass={'build_ext': build_ext},
      zip_safe=False)
