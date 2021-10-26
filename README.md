# cgm-brush
Baryonic post-processing of N-body simulations


Tested on python 3.7.1 and 3.8.5

## Installation Instructions
CGMbrush is not yet available as a public package. For now, you can clone this repository and follow the developer instructions below.

When publically available, you may run the following to install cgmbrush:
$ pip install cgmbrush

## Developer Instructions
If you are forking this repo and intend on working directly with its sources, the following workflow is recommended:

First, setup a virtual environment in the root folder of the respository. Open a terminal and navigate to the base folder of the respository. Type these commands:
$ virtualenv dev_env
$ source dev_env/bin/activate
$ pip install -e .
$ make
$ ipython kernel install --user --name='cgmbrush-dev'

This will setup a virtual environment in the root folder of the respository and install cgmbrush and its requirements in this venv.

By default, output files will be created in '~/cgmbrush/var' and CGMBrush will look for input files (simulation halo tables, etc) from '~/cgmbrush/sims'. Unless this is desirable, you should set environment variables for CGMB_VAR_DIR and CGMB_SIMS_DIR. This could be done by adding "export CGMB_VAR_DIR='...'" to .bashrc, for example.

Then, open tests.ipynb with the ipython kernal you just setup and run the file. All the tests should pass.