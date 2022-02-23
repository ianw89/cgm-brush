###################################################################################################
#
# constants.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################
import numpy as np

kpc = 3.086*10**21 # cm
Mpc = 3.086*10**24 # cm
msun =1.9891 * 10**30 #kilograms
mprot = 1.6726219 * 10**-27 #kilograms


lightspeed = 3e5 #km/s
pi = np.pi
mp = 1.6726216e-24  #proton mass in gr

Yhe =.25 #Helium mass fraction
nPS = msun/mprot*(1-Yhe/2) # accounts for the fact that some mass is in helium 
mu = 1/((1-Yhe) + 2.*Yhe/4)  #mean molecular weight to give number of electrons assuming helium doubly ionized

RadtoSec= 206265 #radians to arsec
MpcInPc = 1e6 
