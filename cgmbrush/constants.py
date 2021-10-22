###################################################################################################
#
# constants.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################
import numpy as np

Mpc = 3.096*10**24 # cm
msun =1.9891 * 10**30 #kilograms
mprot = 1.6726219 * 10**-27 #kilograms
Yhe =.25
nPS = msun/mprot*(1-Yhe/2) # accounts for the fact that some mass is in helium
lightspeed = 3e5 #km/s
pi = np.pi
mp = 1.6726216e-24  #proton mass in gr
            
KPCTOCM = 3.086e21
MPCTOCM = 3.086e24
mu = 1/((1-Yhe) + 2.*Yhe/4)  #mean molecular weight to give number of electrons assuming helium doubly ionized

