# cgm-brush

Baryonic post-processing of N-body simulations. 

Presently CGMBrush is hardcoded to work with data from the Bolshoi Simulations. In the near future it will be organized so that you can easily work with data from other N-body simulations.

Tested on various versions of Python 3. Not tested on Python 2.

## Installation Instructions

CGMbrush is not yet available as a public package. For now, you can clone this repository and follow the developer instructions below.

When publically available, you may run the following to install cgmbrush:
```
$ pip install cgmbrush
```
## Developer Instructions
If you are forking this repo and intend on working directly with its sources, setup instruction are here. You want to install CGMBrush in editable mode (pip install -e path_to_cgmbrush) so that changes made to the sources are picked by files that import it. Because of this, you may want to use an virtual environment so later on if you install a non-developer CGMBrush to your system's python you do not run into issues.

### Setup without a virtual environment, using system python

If you want to setup cgmbrush with the default python on your system, then you only need to install CGMBrush's requirements and CGMBrush itself in editable mode (-e). Open a terminal and navigate to the base folder of the respository. Type these commands:
```
$ pip install -e .
$ pip install -r requirements.txt
```

### Virtual Environment Setup
If you know what a virtual environement is and want to use one, you can follow these instructions to setup a virtual environment in the root folder of the respository. Open a terminal and navigate to the base folder of the respository. Type these commands:
```
$ python -m venv dev_env
$ source dev_env/bin/activate
$ pip install -e .
$ pip install -r requirements.txt
$ ipython kernel install --user --name='cgmbrush-dev'
```
This will setup a virtual environment in the root folder of the respository and install cgmbrush and its requirements in this venv. It also sets up a dedicated ipython kernal for jupyter.



## Input and Output Files

By default, output files will be created in "\~/cgmbrush/var" and CGMBrush will look for input files (simulation halo tables, etc) from "\~/cgmbrush/sims". Unless this is desirable, you should set environment variables for CGMB_VAR_DIR and CGMB_SIMS_DIR. On linux or mac this could be done by adding "export CGMB_VAR_DIR='...'" to .bashrc, and on Windows via the System Properties / Environment Variables... dialog. Alternatively, you could update the values in settings.py of the cgmbrush library.


## Checking your work

Presently, cgmbrush has a tests.ipynb notebook that can be used to gain some confidence that cgmbrush is performing correctly. Currently these tests utilize data from the Bolshoi simulations. If you want to be able to run all of these tests or use cgmbrush with this data, you will need to visit https://www.cosmosim.org/cms/simulations/bolshoi/ and download the Dens256_z0 and BDMV (for z=0) tables. These files should be saved off to your sims directory (see 'Input and Output Files' above) as 'dens256-z-0.csv.gz' and 'halo-z-0.csv.gz'

When you think all is done, open tests.ipynb with the ipython kernel you setup (or your system one that you installed to) and run the file. All the tests should pass.
