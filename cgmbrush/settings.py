###################################################################################################
#
# settings.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################
import os

VAR_DIR = "~/cgmbrush/var" 
"""The directory in which CGMBrush stores files by default, such as the results from 
the various functions of CGMBrush. It is advisable to set this to a location with ample storage.
This can be done with by setting the environment variable CMB_VAR_DIR or by setting the VAR_DIR 
global variable directly in python."""

var_from_env = os.getenv('CGMB_VAR_DIR')
if var_from_env is not None:
    VAR_DIR = var_from_env

SIMS_DIR= "~/cgmbrush/sims" 
"""The directory in which CGMBrush will look for simulation data (density fields, halo tables, etc).
May be set by setting the environment variable CMB_SIMS_DIR or by setting the SIMS_DIR 
global variable directly in python."""

sims_from_env = os.getenv('CGMB_SIMS_DIR')
if sims_from_env is not None:
    SIMS_DIR = sims_from_env
