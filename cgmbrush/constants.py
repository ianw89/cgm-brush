###################################################################################################
#
# constants.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################
import numpy as np

kpc = 3.086*10**21
"""How many centimeters in a Kiloparsec."""
Mpc = 3.086*10**24
"""How many centimeters in a Megaparsec."""
msun =1.9891 * 10**33 
"""How many grams in a Solar mass."""       
mprot = 1.6726219 * 10**-24
"""Proton mass in grams."""       
lightspeed = 2.99792e5
"""Speed of light in a vacuum in kilometers / second."""       

Yhe =.25 
"""Helium mass fraction."""
nPS = msun/mprot*(1-Yhe/2) 
""" accounts for the fact that some mass is in helium """
mean_molecular_weight_electrons = mu = 1/((1-Yhe) + 2.*Yhe/4)  
"""mean molecular weight to give number of electrons assuming helium doubly ionized"""
mean_molecular_weight_baryons = 1/((1-Yhe) + Yhe/4)  
"""the mean molecular weight of just baryons, no"""
fhydrogen = (1-Yhe)/((1-Yhe) + Yhe/4)  
"""fraction of atoms that are hyrogen"""

RadtoSec= 206265 
"""Radians to arsec."""
PcinMpc = 1e6 
"""How many parsecs in a Megaparsec."""
KPCINMPC = 1e3
"""How many kiloparsecs in a Megaparsec."""

DEFAULT_MIN_MASS = 10**10
DEFAULT_MAX_MASS = 8.3*10**14
DEFAULT_MASS_BIN_COUNT = 60 
