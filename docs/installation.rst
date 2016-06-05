Installation Instructions
=========================

.. warning::

  This package has a lingering dependency on an IRAF task that is accessed
  through PyRAF.  IRAF/PyRAF can be obtained through
  `Astroconda <http://astroconda.readthedocs.io/en/latest/>`_.  Though many
  tasks will work without IRAF installed, some will fail under certain
  circumstances.  This depencency will be removed in a future release as
  resources permit.

Install via Anaconda
--------------------

.. note::
  If you do not have Anaconda, please follow the `instructions here
  <https://www.continuum.io/downloads>`_ to install it, or scroll down for
  manual installation of refstis.

After you have anaconda setup, then you can install refstis by
specifying the channel in your install command::

    $ conda install --channel justincely refstis

If you do not yet have IRAF or PyRAF installed, you can install them via the
STScI supplied `AstroConda channel. <http://astroconda.readthedocs.io/>`_

Install from source
-------------------

Refstis can also be installed manually using the source code::

    $ git clone https://github.com/spacetelescope/refstis.git
    $ cd refstis
    $ python setup.py install
