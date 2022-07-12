###################################################################################################
#
# cosmology.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################
import numpy as np
import scipy.integrate as integrate
from cgmbrush.constants import *

class cosmology():
    """
    Cosmological parameters
    """
    h = .7
    OmegaM = 0.27
    OmegaB = 0.0469
    OmegaL = 1- OmegaM #assumes flat universe
    rho_c = 9.31000324385361 *10**(-30) * (Mpc**3 / msun) * (h/.7)**2 # Msun/Mpc^3 from #pcrit = 9.31000324385361e-30 g / cm3
    rho_m = OmegaM*rho_c #density in matter
    #fd = 1 # fraction of baryons in diffuse ionized gas (FRB paper)
    fb = OmegaB/OmegaM #fraction of matter in baryons
    
    @classmethod
    def rhoB(self, z):
        return self.OmegaB*self.rho_c*(1+z)**3

    @classmethod
    def rhoM(self, z):
        return self.OmegaM*self.rho_c*(1+z)**3

    @classmethod
    def elecD(self, z):
        """ 
        Number of electrons per cm^3 assuming helium contributes two electrons.
        """
        return self.rhoB(z)*(1-Yhe/2)/(mprot/msun)/Mpc**3



    ########################################
    # Hubble expansion and distances
    ########################################

    #Hubble function in flat LCDM
    @classmethod
    def Hz(self, z):
        return 100*self.h*(self.OmegaM*(1+z)**3 + self.OmegaL)**.5

    #conformal distance or comoving angular diameter distance
    @classmethod
    def Dconf(self, z):
        return lightspeed*integrate.romberg(self.dConfDistdz,0,z) #integrates a function with given variables

    @classmethod
    def dConfDistdz(self, z):
        return 1./(self.Hz(z))

    #angular diameter distance
    @classmethod
    def DA(self, z):
        return self.Dconf(z)/(1+z)

    #angular size
    @classmethod
    def AS(self, z, physical_radius):
        return physical_radius/(self.DA(z))

    @classmethod
    def comoving_horizon(self, z):
        return lightspeed*integrate.quad(self.dConfDistdz,z,np.inf)[0]

    
class halo():
    """
    Properties of halos.
    """

    #def __init__(self, cosmology):
    #    self.cosmo = cosmology
    #    # z = 0
    #    self.whatthisfor = 0  #

    # 3D NFW profile with r
    # epsilon is a softening
    @classmethod
    def NFW_profile(self, r,rho_nought,R_s, epsilon=0):
        return rho_nought/(((epsilon+r)/R_s)*((1+(r/R_s))**2))    
    
    # Concentration c of a halo from fit to simulation #Matt: we need to cite reference!!!!
    @classmethod
    def halo_conc(self, cosmo, redshift,halo_mass):
        Mpivot = 5.*10**12
        scaling = -0.13
        return (9./(1+redshift))*(halo_mass/Mpivot)**scaling 

    # function for rho_0 of NFW
    @classmethod
    def rho_0(self, cosmo, redshift,halo_mass,R_s):
        c=self.halo_conc(cosmo, redshift,halo_mass)
        return halo_mass/ (4*np.pi*R_s*(np.log(1+c) - c/(1+c)))

    # used to compute virial radius of a halo (important when dark energy is present)
    @classmethod
    def q(self, cosmo, z):
        return cosmo.OmegaL / ((cosmo.OmegaM*(1+z)**3)+ cosmo.OmegaL)  

    # Virial Radius is the critical density of the universe at the given redshift times an 
    # overdensity constant Delta_c using Bryan and Norman (1998), eqn 6, prescription
    @classmethod
    def rho_vir(self, cosmo, z):
        return (18.*np.pi**2 - 82.*self.q(cosmo, z) - 39.*self.q(cosmo, z)**2)*(cosmo.rho_c*(cosmo.OmegaL + cosmo.OmegaM *(1+z)**3))


    #virial density of a halo in solar mass/(physical Mpc)^3
    @classmethod
    def Rvir_den(self,cosmo, z):
        return (4./3. * np.pi * self.rho_vir(cosmo, z))**(1./3.) 


    @classmethod
    def comoving_rvir(self, cosmo, Mhalo, z):
        """Get the co-moving radius for a halo of a given mass at a given redshift usign the Bryan+Norman '98 definition."""
        return (1+z) * (np.power(Mhalo, 1/3) / self.Rvir_den(cosmo, z))

    #comoving radius that is 200 times the matter density in Mpc  (very similar to comoving_rvir; not used by code)
    @classmethod
    def r200Mz(self, cosmo, Mhalo, z):
        rhomatter = cosmo.rho_c*cosmo.OmegaM
        return (Mhalo/((4*np.pi/3)*200.*rhomatter))**(1./3.)
