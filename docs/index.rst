Welcome to radiosonde-teamxuk's documentation!
==============================================

**radiosonde-teamxuk** is a Python package for processing radiosonde data from the TEAMxUK campaign. It converts EDT format files from Vaisala RS41-SGP radiosondes to NCAS-AMOF NetCDF format and generates quicklook plots for data quality assessment.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   usage
   api
   contributing

Features
--------

* Convert EDT format radiosonde data to NCAS-AMOF NetCDF format
* Generate comprehensive quicklook plots
* Calculate derived meteorological variables
* Optional stability analysis plots
* Command-line tools for easy processing

Installation
------------

Quick install from GitHub::

   pip install git+https://github.com/longlostjames/radiosonde-teamxuk.git

For development::

   git clone https://github.com/longlostjames/radiosonde-teamxuk.git
   cd radiosonde-teamxuk
   pip install -e .

See :doc:`installation` for detailed installation instructions.

Quick Start
-----------

Process radiosonde data::

   process-teamxuk-radiosondes input_dir/ output_dir/

Generate quicklook plots::

   generate-teamxuk-radiosonde-quicklooks netcdf_dir/ --stability

See :doc:`quickstart` for more examples.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
