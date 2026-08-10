Usage
=====

Command-line
------------

Many individual tasks have been pulled out into individual modules which
are copied into your path when the package install is run. This allows
a particular procedure to be run on a collection of data without going
through the entire pipeline, and without the limitations that it imposes
(darks/biases from specific proposals, specifically split months/weeks, etc).

Individual tasks for creating super-biases and super-darks may be found below.
These are accessible through the command line.  For now, reference files
delivered to CRDS are now created via an external package, ``run_refstis``,
which handles the logic of which files to combine.

Basejoint
~~~~~~~~~

Basejoint is a confusingly named script that actually creates a baseline bias
file for use in creating baseline dark files.

.. code-block:: bash

  basejoint *.fits outname.fits


Refbias
~~~~~~~

Refbias is the preferred way to create the weekly bias file.  In some cases,
weekbias may be run instead.

.. code-block:: bash

  refbias *.fits outname.fits


Weekbias
~~~~~~~~

If the number of datasets is fewer than a set threshold the weekbias procedure is
run instead of the refbias.  A main difference between the two is that weekbias
uses the baseline bias in some of the calculations.

.. code-block:: bash

  weekbias *.fits outname.fits basebias.fits


Basedark
~~~~~~~~

The basedark task is designed to create a monthly baseline dark.  This is a
combination of all darks for a month that will be used as input when creating
a weekly dark by the weekdark task.

.. code-block:: bash

  basedark *raw.fits outname.fits basebias.fits


Weekdark
~~~~~~~~

The weekdark uses ~a week's worth of data, along with the monthly dark produced
by basedark, to create a dark appriate to a particular week's worth of data.:

.. code-block:: bash

  weekdark *raw.fits outname.fits basedark.fits basebias.fits
