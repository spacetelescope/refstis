Usage
=====

comming soon to documentation near you!

Command-line
------------

Scripting
---------

The Pipeline
------------

The configure file
~~~~~~~~~~~~~~~~~~

The pipeline functionality of the refstis package needs to know some things
about you and where to put stuff.  This is accomplished by parsing a config
file that is assumed to live at ~/config.yaml.

The necessary contents of the file are shown below, though the actuall content
has been omitted.

.. code-block:: yaml

  #  Directories to read/write
  products_directory :
  retrieve_directory :
  delivery_directory :

  # config for querying MAST for data
  mast_server :
  mast_database :
  mast_account :
  mast_password :
  dads_host :

  # config for retrieving from MAST
  archive :
  archive_user :
  email :
  ftp_user :
  host :
