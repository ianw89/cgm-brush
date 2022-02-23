import numpy as np
import scipy.integrate as integrate
from constants import *
# Cosmological parameters




################################################################

###############################################
class cosmology():
    h = .7
    rho_c = 9.31000324385361 *10**(-30) * (Mpc**3 / (msun*1000))*(h/.7)**2 # Msun/Mpc^3 from #pcrit = 9.31000324385361e-30 g / cm3
    OmegaM = 0.27
    OmegaB = 0.0469
    OmegaL = 1- OmegaM #assumes flat universe
    rho_m = OmegaM*rho_c #density in matter
    #fd = 1 # fraction of baryons in diffuse ionized gas (FRB paper)
    fb = OmegaB/OmegaM #fraction of matter in baryons
    
   
    def __init__(self):
        # z = 0
        self.h = h
        
    ########################################
    # Parameters, Constants, and Functions
    ########################################
    def rhoB(self, z):
        return OmegaB*rho_c*(1+z)**3

    def rhoM(self, z):
        return OmegaM*rho_c*(1+z)**3

    # number of electrons per/cm^3 assuming helium contributes two electrons
    def elecD(self, z):
        return fd*rhoB(z)*(1-Yhe/2)/(mprot/msun)/Mpc**3



    ########################################
    # Hubble expansion and distances
    ########################################

    #Hubble function in flat LCDM
    def Hz(self, z):
        return 100*self.h*(self.OmegaM*(1+z)**3 + self.OmegaL)**.5

    #conformal distance or comoving angular diameter distance
    def Dconf(self, z):
        return lightspeed*integrate.romberg(self.dConfDistdz,0,z) #integrates a function with given variables


    def dConfDistdz(z):
        return 1./(Hz(z))

    #angular diameter distance
    def DA(z):
        return Dconf(z)/(1+z)

    #angular size
    def AS(z, physical_radius):
        return physical_radius/(DA(z))

    def comoving_horizon(z):
        return lightspeed*integrate.quad(dConfDistdz,z,np.inf)[0]

from cosmology import cosmology as cosmo

########################################
# Properties of halos
########################################
class halo(cosmo):
    
    # Concentration c of a halo from fit to simulation #Matt: we need to cite reference!!!!
    def halo_conc(self, redshift,halo_mass):
        Mpivot = 5.*10**12
        scaling = -0.13
        return (9./(1+redshift))*(halo_mass/Mpivot)**scaling 

    # function for rho_0 of NFW
    def rho_0(self, redshift,halo_mass,R_s):
        c=halo_conc(redshift,halo_mass)
        return halo_mass/ (4*self.pi*R_s*(np.log(1+c) - c/(1+c)))

    # used to compute virial radius of a hal (important when dark energy is present)
    def q(self, z):
        return cosmo.OmegaL/ ((cosmo.OmegaM*(1+z)**3)+ cosmo.OmegaL)  

    # Virial Radius is the critical density of the universe at the given redshift times an 
    # overdensity constant Delta_c using Bryan and Norman (MM: DATE) prescription
    def rho_vir(self, z):
        return (18.*np.pi**2 - 82.*self.q(z) - 39.*self.q(z)**2)*(cosmo.rho_c*(cosmo.OmegaL + cosmo.OmegaM *(1+z)**3))


    #virial density of a halo in solar mass/(physical Mpc)^3
    def Rvir_den(z):
        return (4./3. * np.pi * self.rho_vir(z))**(1./3.) 

    def comoving_radius_for_halo(Mhalo, z):
        """Get the co-moving radius for a halo of a given mass at a given redshift usign the Bryan+Norman '98 definition."""
        return (1+z)*((Mhalo)**(1./3.) / Rvir_den(z))

    #radius that is 200 times the matter density in Mpc  (very similar to rvir; not used by code)
    def r200Mz(Mhalo, z):
        rhomatter = cosmo.rhocrit*cosmo.OmegaM
        return (Mhalo/((4*np.pi/3)*200.*rhomatter))**(1./3.)
