Installation Instructions
=========================

.. Install via Anaconda
   --------------------
   
   .. note::
     If you do not have Anaconda, please follow the `instructions here
     <https://www.continuum.io/downloads>`_ to install it, or scroll down for
     manual installation of `refstis`.
   
   After you have anaconda setup, then you can install refstis by
   specifying the channel in your install command::
   
       $ conda install --channel XXXXXXXXX refstis
   
   If you do not yet have IRAF or PyRAF installed, you can install them via the
   STScI supplied `AstroConda channel. <http://astroconda.readthedocs.io/>`_

Install from source
-------------------

`Refstis` can be installed manually using the source code::

    $ git clone https://github.com/spacetelescope/refstis.git
    $ cd refstis
    $ python setup.py install
