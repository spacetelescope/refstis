Usage
=====

Command-line
------------

Many individual tasks have been pulled out into individual modules which
are copied into your path when the package install is run. This allows
a particular procedure to be run on a collection of data without going
through the entire pipeline, and without the limitations that it imposes
(darks/biases from specific proposals, specifically split months/weeks, etc).

Basejoint
~~~~~~~~~

Basejoint is a confusingly named script that actually creates a baseline bias
file for use in creating baseline dark files.

.. code-block:: bash

  $ basejoint *.fits outname.fits


Refbias
~~~~~~~

Refbias is the preferred way to create the weekly bias file.  In some cases,
weekbias may be run instead.

.. code-block:: bash

  $ refbias *.fits outname.fits


Weekbias
~~~~~~~~

If the number of datasets is fewer than a set threshold the weekbias procedure is
run instead of the refbias.  A main difference between the two is that weekbias
uses the baseline bias in some of the calculations.

.. code-block:: bash

  $ weekbias *.fits outname.fits basebias.fits


Basedark
~~~~~~~~

The basedark task is designed to create a monthly baseline dark.  This is a
combination of all darks for a month that will be used as input when creating
a weekly dark by the weekdark task.

.. code-block:: bash

  $ basedark *raw.fits outname.fits basebias.fits


Weekdark
~~~~~~~~

The weekdark uses ~a week's worth of data, along with the monthly dark produced
by basedark, to create a dark appriate to a particular week's worth of data.:

.. code-block:: bash

  $ weekdark *raw.fits outname.fits basedark.fits basebias.fits


The Pipeline
------------

.. warning::

  This section is deprecated.  A new set of scripts to automatically derive pipeline 
  reference files is being developed and will be integrated here soon.

The pipeline can be run from any directory, with input/output directories
determined from the config file outlined in the section below.  Simply
executing the pipeline command from the terminal will kick everything off:

.. code-block:: bash

  $ refstis_pipeline

This will perform all the necessary pipeline steps to create the full suite
of superdarks and superbiases, includng;

1. Checking for all the STIS anneal datasets and determining anneal months.
2. Retrieving new raw dark and bias observations for each month.
3. Reducing the raw datasets and combining into weekly and bi-weekly reference files.
4. Prepping the reference files for delivery into CRDS.
5. Running the data against a test-suite and verifing no errors are produced.

After all steps have been completed for a given anneal month, the pipeline will
send you an email with the delivery form.  This should be forwarded to
redcat@stsci.edu once any final checks on the data are perfomed.

The configure file
~~~~~~~~~~~~~~~~~~

The pipeline functionality of the refstis package needs to know some things
about you and where to put stuff.  This is accomplished by parsing a config
file that is assumed to live at ~/refstis_config.yaml.

The necessary contents of the file are shown below, though the content is
dummy and will need to be configured for you specifically.

.. code-block:: yaml

  #  Directories to read/write
  products_directory : '/Users/myself/refstis/data/'
  retrieve_directory : '/Users/myself/refstis/requested/'
  delivery_directory : '/Users/myself/refstis/to_deliver/'


  # config for querying MAST for data
  mast_server : 'server@name.stsci.edu'
  mast_database : 'db_name'
  mast_account : 'username'
  mast_password : 'Pa$$w3rD'
  dads_host : 'dads_host.stsci.edu'

  # config for retrieving from MAST
  archive : 'archive.stsci.edu'
  archive_user : 'myself'
  email : 'myself@stsci.edu'
  ftp_user : 'myself'
  host : 'host.domain.com'

  # Proposals to use for darks/biases
  dark_proposals:
   - 12000
   - 13001
   - 14243

  bias_proposals:
   - 12001
   - 13005
   - 14244
