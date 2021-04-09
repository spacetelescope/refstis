#!/usr/bin/env python

from setuptools import setup, find_packages
import os
import glob

setup(
    name = 'refstis',
    version = '0.8.1a1',
    description = 'Pipeline to create STIS CCD superdarks and superbiases',
    author = 'Allyssa Riley, Sean Lockwood, Justin Ely',
    author_email = 'https://hsthelp.stsci.edu',
    url = 'https://refstis.readthedocs.io',
    packages = find_packages(),
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python :: 3',
                   'Development Status :: 4 - Beta',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    scripts = glob.glob('scripts/*'),
    install_requires = ['pyyaml',
                        'numpy',
                        'astropy>=3.1',
                        'stistools'],
    )
