#!/usr/bin/env python

from setuptools import setup, find_packages
import os
import glob

setup(
    name = 'refstis',
    version = '0.6.1',
    description = 'Pipeline to create STIS CCD superdarks and superbiases',
    author = 'Justin Ely',
    author_email = 'ely@stsci.edu',
    packages = find_packages(),
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    scripts = glob.glob('scripts/*'),
    install_requires = ['pyyaml',
                        'numpy>=1.10',
                        'astropy>=1.0.1',
                        'stistools>=1.0.2'],
    )
