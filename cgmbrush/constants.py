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
mean_molecular_weight_electrons = mu = 1/((1-Yhe) + 2.*Yhe/4)  #mean molecular weight to give number of electrons assuming helium doubly ionized
mean_molecular_weight_baryons = 1/((1-Yhe) + Yhe/4)  #the mean molecular weight of just baryons, no
fhydrogen = (1-Yhe)/((1-Yhe) + Yhe/4)  #fraction of atoms that are hyrogen

#print(mu, mean_molecular_weight_baryons, fhydrogen)

RadtoSec= 206265 #radians to arsec
MpcInPc = 1e6 
