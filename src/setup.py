#!/usr/bin/env python
import os

"""
setup.py file for proxTV python interface
"""

from distutils.core import setup, Extension


proxtv_module = Extension('_proxtv_internal',
                           sources=['proxtv_wrap.cxx', 'proxtv.cc'],
                           library_dirs=[os.getcwd()],
                           libraries=['proxtv','blas','lapack'],
                           extra_compile_args=['-fopenmp'],
                           extra_link_args=['-lgomp']
                           )

setup (name = 'proxtv_internal',
       version = '0.1',
       author      = "Alvaro Barbero, Suvrit Sra",
       description = """Python interface por proxTV""",
       ext_modules = [proxtv_module],
       py_modules = ["proxtv_internal"],
       )

