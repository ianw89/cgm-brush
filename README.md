# CGMBrush: what is it?

CGM Brush is a python library for baryonic post-processing of N-body simulations. CGM Brush was initially built with the goal of redistributing mass from dark matter only cosmological simulations in accordance with various models of the circumgalactic medium (CGM). The routines work in 2D and utiltze fast convolutions, which allows the library to work very quickly, and allows the the library to work at much higher resolutions than the original cosmological simulations grid.

This tool can be used for a variety of purposes, but the initial application is using fast radio bursts to study and constrain the CGM.

Presently CGMBrush is hardcoded to work with data from the Bolshoi Simulations. In the near future it will be organized so that you can easily work with data from other N-body simulations.

CGMBrush is Tested on various versions of Python 3 and sci-py the latest scipy packages as of 2021. It is not tested on Python 2.

## Installation Instructions

CGMbrush is not yet available as a public package. For now, you can clone this repository and follow the developer instructions below.


## Developer Instructions
If you are forking this repo and intend on working directly with its sources, setup instruction are here. You want to install CGMBrush in editable mode (pip install -e path_to_cgmbrush) so that changes made to the sources are picked by files that import it. Because of this, you may want to use an virtual environment so later on if you install a non-developer CGMBrush to your system's python you do not run into issues.

### Setup without a virtual environment, using system python

If you want to setup cgmbrush with the default python on your system, then you only need to install CGMBrush's requirements and CGMBrush itself in editable mode (-e). Open a terminal and navigate to the base folder of the respository. Type these commands:
```
$ pip install -e .
$ pip install -r requirements.txt
```

### Virtual Environment Setup
If you know what a virtual environement is and want to use one, you can follow these instructions to setup a virtual environment in the root folder of the respository instead. Open a terminal and navigate to the base folder of the respository. Type these commands:
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


## Bolshoi Simulations Files
CGM Brush uses .csv.gz density field and halo table files from the Bolshoi simulations. You will need to visit https://www.cosmosim.org/cms/simulations/bolshoi/ and read and understand the Dens and BDMProf tables. The density field files are acquired by querying for Dens256_z0, Dens512, and/or Dens256 (for whatever redshifts). For instance, you go to https://www.cosmosim.org/query and enter
```
SELECT * FROM Bolshoi.Dens256_z0
```
to get the density field at 256x256 resolution for z=0. You also need the halo tables (BDMProf), which should be acquired and saved off seperatly for each redshift snapshot desired. After you queries complete, you can download a .csv file from the cosmosim job and then run gzip on it. These files should be saved off and renamed your sims directory (see 'Input and Output Files' above) as 'dens256-z-0.0.csv.gz' and 'halo-z-0.0.csv.gz', for instance. Naming of higher reshift files is more complicated; see the comments in the BolshoiProvider class in the cgmbrush.py file.


## Checking your work

Presently, cgmbrush has a tests.ipynb notebook that can be used to gain some confidence that cgmbrush is performing correctly. Currently these tests utilize data from the Bolshoi simulations, so you must have some of these files available as described above.

When you think you are ready, open tests.ipynb with the ipython kernel you setup (or your system one that you installed to) and run the file. If you've downloaded and set up the z=0 density field and halo table correctly, most of the tests will pass (a few require higher redshift files, and this is noted in the tests comments). After you are satisfied, we suggest that you then proceed to the tutorial.ipynb to learn a little more about how to use CGM Brush.

## Use and citations
If you are using CGM Brush for scientific applications, we encourage you to contact ianw89@live.com, adnank@uw.edu, or mcquinn@uw.edu.

The paper detailing the methods used in this library is not yet published under a peer-reviewed journal. Please check back later for details on how to cite your use of CGM Brush.
