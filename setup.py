from distutils.core import setup
import os
import glob

setup(
    name = 'refstis',
    version = '0.0.1', 
    description = 'Replacement pipeline for STIS superdarks and superbiases',
    author = 'Justin Ely',
    author_email = 'ely@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    scripts = glob.glob('scripts/*'),
    requires = ['numpy','pyfits'],
    )
