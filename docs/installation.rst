Installation Instructions
=========================

.. warning::

  This package has a lingering dependency on an IRAF task that is accessed
  through PyRAF.  IRAF/PyRAF can be obtained through
  `Astroconda <http://astroconda.readthedocs.io/en/latest/>`_.  Though many
  tasks will work without IRAF installed, some will fail under certain
  circumstances.  This depencency will be removed in a future release as
  resources permit.

Install from source
-------------------

Refstis can also be installed manually using the source code::

    $ git clone https://github.com/spacetelescope/refstis.git
    $ cd refstis
    $ python setup.py install
