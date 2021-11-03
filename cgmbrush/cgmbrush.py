###################################################################################################
#
# cgmbrush.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################

from __future__ import print_function 
from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.integrate as integrate
import random
from numpy.polynomial import Polynomial as P
import pandas as pd
import os
import abc
import cProfile
import io
import pstats
import datetime
import math
import sys
from scipy.interpolate import interp1d

from cgmbrush.settings import *
from cgmbrush.constants import *

########################################
# Parameters, Constants, and Functions
########################################

# Cosmological parameters
z = 0
h = .7
rho_c = 9.31000324385361 *10**(-30) * (Mpc**3 / (msun*1000))*(h/.7)**2 # M./Mpc^3 from #pcrit = 9.31000324385361e-30 g / cm3
OmegaM = 0.27
OmegaL = 1- OmegaM #assumes flat universe
OmegaB = 0.0469 
rho_m = OmegaM*rho_c
fd = 1 # fraction of baryons in diffused ionized gas (FRB paper)
fb = OmegaB/OmegaM #fraction of matter in baryons

def rhoB(z):
    return OmegaB*rho_c*(1+z)**3

def rhoM(z):
    return OmegaM*rho_c*(1+z)**3

# units: #/cm^3
def elecD(z):
    return fd*rhoB(z)*(1-Yhe/2)/(mprot/msun)/Mpc**3




########################################
# Configuration of the simulation
# ########################################

# TODO Bolshoi specific; need to make these properties of the SimulationProvider
L=250/h # length of box in Mpc
Ncel= 256  # number of cells of density field
dx= L/Ncel # width of each cell
V= dx**3 # volume of each cell



########################################
# Properties of halos
########################################

# Concentration c of a halo:
def halo_conc(redshift,halo_mass):
    return (9/(1+redshift))*(halo_mass/(5*10**12))**-0.13  #MM added the exponent

# function for rho_0 of NFW
def rho_0(redshift,halo_mass,R_s):
    c=halo_conc(redshift,halo_mass)
    return halo_mass/ (4*pi*R_s*(np.log(1+c) - c/(1+c)))

# virial radius of a halo: 
def q(z):
    return OmegaL/ ((OmegaM*(1+z)**3)+ OmegaL)  

# Virial Radius is the critical density of the universe at the given redshift times an 
# overdensity constant Δ_c. Several Δ_c conventions exist; some are redshift dependent.
def rho_vir(z):
    return (18*pi**2 + 82*q(z) - 39*q(z)**2)*(rho_c*(OmegaL + OmegaM *(1+z)**3))

def Rvir_den(z):
    return (4/3 * np.pi * rho_vir(z))**(1/3) # physical, units in 1/r**3

def comoving_radius_for_halo(Mhalo, z):
    """Get the co-moving radius for a halo of a given mass at a given redshift usign the Bryan+Norman '98 definition."""
    return (1+z)*((Mhalo)**(1/3) / Rvir_den(z))

#radius that is 200 times the matter density in Mpc  (very similar to rvir)
def r200Mz(Mhalo, z):
    rhocrit = 2.7755e11*h*h
    rhomatter = rhocrit*OmegaM
    return (Mhalo/((4*np.pi/3)*200.*rhomatter))**(1./3.)



########################################
# Histogram functions
########################################

# The function takes an array, bins, and number of cells, and creates a histogram

# Histogram array with min and max defined
def histArray(arr,nbin,Ncel,DMmin,DMmax):
    hist, bins = np.histogram(arr,bins=nbin,range=(DMmin, DMmax))
    dv=bins[1]-bins[0]
    hist=(hist)/(dv*Ncel**2)
    cumulative = (np.cumsum(hist)*dv)
    center = (bins[:-1] + bins[1:]) / 2
    return center, hist, cumulative[-1]


# Histogram array without min and max
# def histArray(arr,nbin,Ncel):
#     hist, bins = np.histogram(arr,bins=nbin)#,range= (np.mean(arr),5*np.mean(arr)))
#     dv=bins[1]-bins[0]
#     hist=(hist)/(dv*Ncel**2)
#     cumulative = (np.cumsum(hist)*dv)
#     center = (bins[:-1] + bins[1:]) / 2
#     return center, hist, cumulative[-1]

# plots the histogram of the dispersion measures
def plothist(arr,nbin,xmin,xmax,title,Ncel):
    hist, bins = np.histogram(arr,bins=nbin)
    dv=bins[1]-bins[0]
    hist=(hist)/(dv*Ncel**2)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.title(title)
    plt.show()

    

def plothist2(arr,nbin,xmin,xmax,title):
    for i in range(len(arr[:,0,0])):    
        hist, bins = np.histogram(arr[i,:,:],bins=nbin)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        #plt.bar(center, hist, align='center', width=width)
        plt.plot(center, hist)
        plt.xlim(xmin=xmin,xmax=xmax)
        plt.title(title)
    plt.legend(loc='upper right')
    plt.show()




########################################
# Function to create arrays with an upper and lower bound for halo mass
########################################

# min, max in solar masses
def haloArray_minmax(pdHalosN, min, max):
    pdHalosN = pdHalosN.drop(pdHalosN[pdHalosN.Mvir < min].index)
    pdHalosN = pdHalosN.drop(pdHalosN[pdHalosN.Mvir > max].index)
    return pdHalosN




########################################
# Cosmology functions
########################################

RadtoSec= 206265; #radians to arsec

def Hz(z):
    return 100*h*(OmegaM*(1+z)**3 + OmegaL)**.5

#conformal distance is comoving distance
def Dconf(z):
    return lightspeed*integrate.romberg(dConfDistdz,0,z) #integrates a function with given variables

def dConfDistdz(z):
    return 1./(Hz(z))

#angular diameter distance
def DA(z):
    return Dconf(z)/(1+z)

#angular size
def AS(z,GRadius):
    return GRadius/(DA(z))


# comoving horizon
def Dconf2(z):
    return lightspeed*integrate.quad(dConfDistdz,z,np.inf)[0]




########################################
# Miscelaneous functions
########################################

# Dispersion measure: analytical function
def DM_analytical(z):
    return lightspeed* integrate.romberg(integrandDM,0,z) #integrates a function with given variables

# integrand of dispersion measure
def integrandDM(z):
    return 10**6*dConfDistdz(z)*elecD(z)/(1+z)**2 

# normalize DM
def normDM(DM,z):
    return DM *10**6* dx* elecD(z) /(1+z)**2  # Psc cm-2

def redshifted_DM(DM,z):
    return DM *(1+z)  # Psc cm-2

# finding redshift for comoving distance
def RSvsCDar(z):
    RSvsCD = np.zeros((100,2))
    ar=np.linspace(0,z,100)
    for i in range(0,len(ar)):
        RSvsCD[i,0] = ar[i]
        RSvsCD[i,1] = Dconf(ar[i])
    return RSvsCD

# function to get redshift from the comoving distance by interpolating upto z, returns Mpc/h 
def CDtoRS(CD,z):
    x = RSvsCDar(z)[:,0]
    y = RSvsCDar(z)[:,1]
    p = P.fit(x, y, 2)
    return ((p - CD).roots()[0])
    
# number of boxes needed for a given redshift
def numBoxes(z):
    return float(round(Dconf(z)/(L)))


# finds electron density for a given box number n of the stack
def elecDforBox(n):
    avgZ = (CDtoRS(n*L,1)+CDtoRS((n-1)*L,1))/2
    return elecD(avgZ)

def avgZ(n): 
    return max((CDtoRS(n*L,1)+CDtoRS((n-1)*L,1))/2,0)

# normalizes DM for a given box and accounts for electron density 
def normDM(DM,z):
    return (DM) *10**6* dx* elecD(z) /(1+z)**2  # Psc cm-2
    
def z_eff(zmin,zmax,L_box):
    return ((DM_analytical(zmax)-DM_analytical(zmin))/((L_box*10**6)*elecD(0)))**(1) - 1

def RS_array_gen(z_max,L):
    """Computes redshift values for each box when stacking boxes of length L out to redshift z_max.
    
    Returns an array of these z values.
    """
    RS_array = []
    RS_array.append(CDtoRS(L,1))
    for i in range(1,int(numBoxes(z_max))):
        RS_array.append(z_eff(CDtoRS(float(L*i),z_max),CDtoRS(float(L*(i+1)),z_max),L))

    return RS_array

    





########################################
# This section is the interface and code for extracting density fields and halos
# from N-body simulations outputs for use in cgmbrush.
########################################

class SimulationProvider(metaclass=abc.ABCMeta):
    """Interface used by cgmbrush to access information from the N-body simulation.
    
    Users are expected to subclass this interface. See examples."""

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_density_field') and
        callable(subclass.get_density_field) and 
        hasattr(subclass, 'get_z_name') and
        callable(subclass.get_z_name) and 
        hasattr(subclass, 'get_halos') and
        callable(subclass.get_halos) or
        NotImplemented)

    @abc.abstractmethod
    def get_density_field(self, redshift: float, resolution: int):
        """Gets the density field for the given redshift and grid resolution."""
        raise NotImplementedError

    @abc.abstractmethod
    def get_halos(self, redshift : int) -> pd.DataFrame:
        """Gets halo information for the given redshift."""
        raise NotImplementedError

    @abc.abstractmethod
    def get_z_name(redshift: float) -> str:
        """Gets the shorthand name associated with a given redshift, used in filenames."""
        raise NotImplementedError


class BolshoiProvider(SimulationProvider):
    """SimulationProvider implementation for Bolshoi simulations.

    Reads density fields for downloaded from the Bolshoi Simulation archives named like dens<resolution-z-<zzz>.csv.gz and halo-z-<zzz>.csv.gz,
    where <resolution> is either 256 or 512 and <zzz> is a shorthand for redshift of the data, which can be related to the actual redshift
    via the z_to_filename map.
    
    Note for dealing with z > 0:
    In order to deal with the fact that there is different density fields and halos for specific redshifts, 
    a filename convention is created here that maps specific redshift values to a shorter filename moniker.
    So, the data for z≈0.31 is in the files with 0.4 in them as per the default z_to_filename map. The map can 
    be re-written as per your application.
    """

    def __init__(self):
        # Associative array of (redshift, resolution) => 3D numpy grid 
        # This works very well with np.savez, which is much faster to read
        # than the Bolshoi density files
        self.density_fields = {}  
        self.halos = {}
        self.z_to_filename = {
            '0':'0',
            str(round(0.0845770240448820, 2)): '0.1', 
            str(round(0.1337050904798766, 2)): '0.2',
            str(round(0.2225162656597088, 2)): '0.3',
            str(round(0.3113821868867406, 2)): '0.4',
            str(round(0.4089786628982832, 2)): '0.5',
            str(round(0.5101330792301793, 2)): '0.6',
            str(round(0.6205902514416375, 2)): '0.7', 
            str(round(0.7454876185015988,2)): '0.8', 
            str(round(0.8931355846491102,2)): '0.9'
        }

    def get_z_name(self, redshift: float) -> str: 
        return self.z_to_filename[str(round(redshift, 2))]

    def get_density_field(self, redshift: float, resolution: int):
        """Gets the density field for the given redshift and grid resolution."""

        if str(round(redshift, 2)) not in self.z_to_filename:
            raise ValueError("There is no file associated with the z={}. Add an entry to the z_to_filname dictionary on this object.".format(redshift))
        z_name = self.get_z_name(redshift)

        # If in memory, use it
        if (z_name, resolution) in self.density_fields:
            return self.density_fields[(z_name, resolution)]

        # If not, try fast reading from saved off numpy array            
        filename = "bol_den_field_" + z_name + '_' + str(resolution) 
        result = []
        try:
            result = loadArray(filename)
        except IOError:
            pass 
        
        # Do the slow reading from the bolshoi file and save off fast copy
        if len(result) == 0:
            result = self.import_density_field(z_name, resolution)
            if len(result) == 0:
                raise IOError("There was a problem importing the Bolshoi density field for (" + str(redshift) + ", " + str(resolution) + ")")
            saveArray(filename, result)

        self.density_fields[(z_name, resolution)] = result
        return result

    def get_halos(self, redshift: float) -> pd.DataFrame:
        """Gets halo information for the given redshift."""

        if str(round(redshift, 2)) not in self.z_to_filename:
            raise ValueError("There is no file associated with the z={}. Add an entry to the z_to_filname dictionary on this object.".format(redshift))
        z_name = self.get_z_name(redshift)

        # If in memory, use it
        if z_name in self.halos:
            return self.halos[z_name]

        # If not, try fast reading from saved off numpy array            
        filename = "bol_halos_" + z_name
        result = []
        try:
            result = loadArray(filename)
            # The rest of the code expects this to be a DataFrame, so convert it back
            result = pd.DataFrame(result, columns = ['row_id','x','y','z','Mvir','Mtot','Rvir','ix','iy','iz'])
        except IOError:
            pass 
        
        # Do the slow reading from the bolshoi file and save off fast copy
        if len(result) == 0:
            result = self.extract_halos(z_name)
            if len(result) == 0:
                raise IOError("There was a problem importing the Bolshoi halos for redshift" + str(redshift))
            saveArray(filename, result)

        self.halos[z_name] = result
        return result

    def import_density_field(self, z_name, resolution):
        """Imports the density field for a given redshift and smooths it. Stores the result in a cache that is faster to access."""
        
        if resolution != 256 and resolution != 512:
            raise ValueError("Only resolution 256 or 512 is supported.")

        resStr = str(resolution)

        # Have to hardcode z=0 table because of the unique column names
        if z_name == '0':
            # reading density field and halos data
            file_path= os.path.join(SIMS_DIR, 'dens'+resStr+'-z-0.csv.gz')
            pdDens=pd.read_csv(file_path)

            # extracting columns
            pdDensN=pdDens[['Bolshoi__Dens'+resStr+'_z0__ix','Bolshoi__Dens'+resStr+'_z0__iy','Bolshoi__Dens'+resStr+'_z0__iz','Bolshoi__Dens'+resStr+'_z0__dens']]

            # 3D density array
            pdDensN=pdDensN.sort_values(['Bolshoi__Dens'+resStr+'_z0__ix','Bolshoi__Dens'+resStr+'_z0__iy','Bolshoi__Dens'+resStr+'_z0__iz'])
            tden = pdDensN['Bolshoi__Dens'+resStr+'_z0__dens'].values
            tden2 = np.reshape(tden,(resolution,resolution,resolution))
            return normDM((tden2+1).sum(2), 0)

        else:
            file_path = os.path.join(SIMS_DIR, 'dens'+resStr+'-z-{}.csv.gz'.format(z_name))
            den = pd.read_csv(file_path)
            den2=den[['Bolshoi__Dens'+resStr+'__ix','Bolshoi__Dens'+resStr+'__iy','Bolshoi__Dens'+resStr+'__iz','Bolshoi__Dens'+resStr+'__dens']]

            # 3D density array
            den_sorted=den2.sort_values(['Bolshoi__Dens'+resStr+'__ix','Bolshoi__Dens'+resStr+'__iy','Bolshoi__Dens'+resStr+'__iz'])
            den_vals = den_sorted['Bolshoi__Dens'+resStr+'__dens'].values
            den = np.reshape(den_vals,(resolution,resolution,resolution))
            
            den = normDM((den+1).sum(2),0) 

            test_l= (np.repeat((np.repeat(den,4,axis=0)),4,axis=1))

            test_sm =  gauss_sinc_smoothing(test_l,4,4,1)
            smoothed_field = test_sm.reshape([256, 4, 256, 4]).mean(3).mean(1)
            return smoothed_field

        
    # Extract halos for a given redshift
    def extract_halos(self, z_name):
        name = 'halo-z-{}.csv.gz'.format(z_name)    
        file_path= os.path.join(SIMS_DIR, name)
        halos = pd.read_csv(file_path)
        return halos











########################################
# Convolution Functions
########################################

def my_convolve(a, b):
    """FFT convolve. Assumes a is the bigger 2D array and b is the smaller 2D mask."""
    halfwidth = int(b.shape[0] / 2)
    a2 = np.roll(np.roll(a, -halfwidth, axis=0), -halfwidth, axis=1)
    return np.fft.irfft2(np.fft.rfft2(a2) * np.fft.rfft2(b, a2.shape))

# Create halo array from halo table for convolution
def create_halo_array_for_convolution(pdHalos, M_min, M_max, logchunks):
    """Creates an array of indexes corresponding to indexes in pdHalos at the edges of logarithmic mass bins.

    The length will be logchunks."""
    halos = haloArray_minmax(pdHalos,M_min,M_max)
    df = halos.sort_values(by='Mvir',ascending=True)
    sorted_haloMasses = df['Mvir'].values
    histH, binsH = np.histogram(sorted_haloMasses,bins=np.logspace(np.log10(M_min),np.log10(M_max), logchunks))
    bins=np.append([0],histH.cumsum()).tolist()

    return df,bins

def project_spherical_3Dto2D(f, x, y, Rvir):
    """Project a spherically symmetric profile from 3D to 2D.
    
    Rvir is in units of virial radii per cell, and is used as a hard boundary of the sphere.""" 
    #print("Rvir, x, y: {}, {}, {}".format(Rvir, x, y))
    boundary = math.sqrt(max(0.0, Rvir**2-(x**2+y**2)))
    if boundary == 0.0:
        return 0.0
    else:
        #print("Boundary: {}".format(boundary))
        return 2 * integrate.quad(f, 0, boundary, args=(x,y))[0] 
    
def NFW2D(x, y, rho_nought, R_s, Rvir):
    """Project NFW profile from 3D to 2D.
    
    Rvir is in units of virial radii per cell.""" 
    offset=float(.1) # TODO .5 and move offset
    boundary = math.sqrt(max(0.0, Rvir**2-(x**2+y**2)))
    if boundary == 0.0:
        return 0.0
    else:
        #print("Boundary: {}".format(boundary))
        return 2 * integrate.quad(lambda x, y, z: rho_nought/(((offset+(x**2+y**2+z**2)**.5)/R_s)*(1+((x**2+y**2+z**2)**.5)/R_s)**2), 0, boundary, args=(x,y))[0]

# Function creates a smaller grid from a larger grid by smoothing
def smoothfield(big, nbig, nsmall):
    
    roll = int((nbig/nsmall)/2)
    y = np.roll(big, roll, axis=0)  #so that y grid takes 2 cells behind and 2 in front
    y = np.roll(y, roll, axis=1)

    small = y.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
    return small

# Function subtracts halo field from density field
def removeConvolvedHalos(densityField,haloField):
    return densityField-haloField

# Gaussian and sinc smoothing of the halo field
def gauss_sinc_smoothing(smoothing_array,sigma_gauss,width_sinc,resolution):
    
    func = smoothing_array
    fft_func = np.fft.fft2(func)

    # sigma is the variance of the gaussian. sigma**2 = 16 or 3*sigma**2 = 16
    sigma = sigma_gauss

    # nsinc is the full width of the box function. So 4 means two boxes on either side of the point.
    nsinc = width_sinc

    # wavenumbers for the transforms
    kx = (np.fft.fftfreq(1024*resolution))
    ky = kx[:,np.newaxis]

    #     ky = (np.fft.fftfreq(1024*resolution))

    prod = np.zeros([1024*resolution,1024*resolution],dtype=complex)

#     for i in range(1024*resolution):
#         for j in range(1024*resolution):

#             # correct functions
#             prod[i,j] = fft_func[i,j]*np.exp(-(sigma**2) *2*np.pi**2* (kx[i]**2 + ky[j]**2 ))
            
#             prod[i,j] =  prod[i,j]*(np.sinc(kx[i]*nsinc*resolution)*np.sinc(ky[j]*nsinc*resolution))**2


    prod = fft_func*np.exp(-(sigma**2) *2*np.pi**2* (kx**2 + ky**2 )) #* (np.sinc(kx*nsinc*resolution)*np.sinc(ky*nsinc*resolution))**2
    result = np.fft.ifft2(prod)
    
    return result.real


class CGMProfile(metaclass=abc.ABCMeta):
    """Interface used by cgmbrush that handles creating a CGM profile to convolve."""

    def __init__(self):
        self.debug = False

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_mask') and
        callable(subclass.get_mask) or 
        NotImplemented)

    @abc.abstractmethod
    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):
        """Constructs and returns a 2D mask to convolve with the halo locations
           for the specified halo parameters and redshift."""
        raise NotImplementedError

class MassDependentProfile(CGMProfile):
    """Profile to use when you want to use different profiles for different mass ranges."""

    def __init__(self, lower_profile: CGMProfile, higher_profile: CGMProfile, cutoff_mass: float):
        self.name = lower_profile.name + "_and_" + higher_profile.name + "_{:.1f}".format(math.log10(cutoff_mass))
        self.lower_profile = lower_profile
        self.higher_profile = higher_profile
        self.cutoff_mass = cutoff_mass
        super().__init__()

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):
        if mass <= self.cutoff_mass:
            return self.lower_profile.get_mask(mass, comoving_radius, redshift, resolution, scaling_radius, cellsize, fine_mask_len)
        else:
            return self.higher_profile.get_mask(mass, comoving_radius, redshift, resolution, scaling_radius, cellsize, fine_mask_len)


class TophatProfile(CGMProfile):

    def __init__(self):
        self.name = "tophat"
        super().__init__()

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):

        # fine mask
        # fine mask size has to correspond to the size of the mask that I eventually trim
        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)
        scale_down = 2  # making the grid coarser    
        
        # TODO Why are we multiplying by scale_down
        r = (x**2+y**2)**.5
        fine_mask = r <= (scaling_radius * scale_down * comoving_radius / cellsize) 
        return fine_mask.astype(float)

class SphericalTophatProfile(CGMProfile):

    def __init__(self, extra=1):
        self.name = "STH"
        self.extra = extra
        super().__init__()

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)
        scale_down = 2  # making the grid coarser    

        r = (x**2+y**2)**.5

        # TODO this 'extra' thing I added was for compatibility with the 2RVSTH_and_NFW_X profiles, which for some reason have an extra *2 in the fine mask line below. 
        # I think it's a BUG, but all the profiles like 2RVSTH_and_NFW_13.5 had it in their spherical tophat code (not the NFW side, curiosuly)
        fine_mask = r <= (self.extra * scaling_radius * scale_down * comoving_radius / cellsize)
        fine_mask=fine_mask.astype(float)
        
        # TODO r should always be less than Rv...
        Rv = (scaling_radius * scale_down * comoving_radius / cellsize)
        fine_mask = fine_mask * ((1-((r/Rv)**2))**(2))**(1/4)
        return fine_mask

class NFWProfile(CGMProfile):
    
    def __init__(self):
        self.name = "NFW"
        super().__init__()

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):

        # TODO ? add redshift to the function above ?
        R_s= comoving_radius/(halo_conc(redshift,mass)*cellsize)  
        rho_nought = rho_0(redshift,mass,R_s)
        scale_down = 2  # making the grid coarser    

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)
        
        if self.debug:
            with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
                print('R_s: {}, rho_0: {}, r/size: {}'.format(R_s, rho_nought, comoving_radius / cellsize))

        vec_integral = np.vectorize(NFW2D)
        fine_mask = vec_integral(x, y, rho_nought, R_s, scale_down * comoving_radius / cellsize)

        if self.debug:
            with np.printoptions(precision=1, linewidth=1000, threshold=sys.maxsize):
                print("from vec integral:")
                print(fine_mask)                
                print("done")

        r=(x**2+y**2)**.5 
        
        fine_mask = fine_mask.astype(float)
        fine_mask[r > scale_down * comoving_radius / cellsize] = 0 # sets very small numbers to 0

        return fine_mask

class FireProfile(CGMProfile):
    """
    This is FIRE simulation profile (from eyeballing https://arxiv.org/pdf/1811.11753.pdf and using their $r^{-2}$ scaling)
    I'm using the $r^{-2}$ profile they find, with an exponetial cutoff at rmax, where we use conservation of mass to determine rmax. 
    
    So: $$\\rho = \\rho_{0}\left(\\frac{r}{r_*} \\right)^{-2} \exp[-r/r_{\\rm max}]$$

    Their results technically hold for $10^{10} - 10^{12} M_\odot$ and $z=0$, but I think it's reasonable to assume they extrapolate to other moderate redshifts and somewhat more massive halos.

    To normalize things we are using that:
    $$M_{\rm gas} = \int 4 \pi r^2 dr \rho_{0} \left(\frac{r}{r_*} \right)^{-2} \exp[-r/r_{\rm max}] = 4 \pi r_{\rm max} \rho_0 r_*^2$$

    Note that sometimes I'm working with number densities and sometimes mass densities.
    """
    
    def __init__(self):
        self.name = "fire"
        super().__init__()

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)             

        Msun =  1.9889e33  # gr 
        adjustmentfactor = 1  #probably we want to range between 0.5-2 as the numbers look sensible in this range
        
        # These are specifications taken from the fire simulations
        RinterpinRvir = 0.3  # this is the point where I read off the density nrmalization
        logMinterp = np.array([10., 11., 12.])  # these are log10 of the halo masses they consider
        nHinterp = np.array([0.5e-4, 0.8e-4, 1e-4])  # these are their number densities in cubic cm
        nHinterp = interp1d(logMinterp, nHinterp, fill_value="extrapolate")

        rho0 = adjustmentfactor * nHinterp(np.log10(mass))
        Rinterp = RinterpinRvir * comoving_radius #rvir(mass, z)
        
        Ntot = mass*Msun*fb/(mu*mp)/(MPCTOCM**3) # TODO msun here is is grams, but elsewhere it is in kg. Double check math.
        rmax = Ntot/(4.*np.pi*rho0*Rinterp**2 )  #from integrating above expression for rho

        # TODO these two lines are dead code, they were returned but not used before. What is the deal?
        #rarr = np.logspace(-2, np.log10(5*rmax), 100) # Cut off at 5 times exponential cutoff
        #rhoarr = rho0*(rarr/Rinterp)**-2*np.exp(-rarr/rmax) #number density: per cm3
        
        ## creating a mask
        # TODO perf bottleneck is here
        # TODO use symmetry to save computation
        f1 = lambda x, y, z: fire_func(((x**2+y**2+z**2)**.5), rmax, Rinterp, rho0, cellsize)
        
        vec_integral = np.vectorize(project_spherical_3Dto2D)   

        fine_mask = vec_integral(f1, x, y, rmax / cellsize) 
        if self.debug:
            with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
                print("From vec integral:")
                print(fine_mask)                

        fine_mask = fine_mask.astype(float)

        return fine_mask
    
# From rhogas Fire, TODO reorganize
# The user can create their own function
# All length scales have to be converted into units of cellsize
def fire_func(r, rmax, Rinterp, rho0, cellsize):
    R1 = rmax/cellsize  # Convert length scale into units of cellsize
    assert r < R1, 'r should always be less than R1, but r={} and R1={}.'.format(r, rmax)
    R2 = Rinterp/cellsize
    #print('r: {}, result: {}'.format(r, rho0 * np.exp(-r/R1) * ((r+.5)/R2)**-2))
    return rho0 * np.exp(-r/R1) * ((r+.5)/R2)**-2

class PrecipitationProfile(CGMProfile):
    """
    Voit Perciptation limited model from Appendix A in https://arxiv.org/pdf/1811.04976.pdf.
    """
    
    def __init__(self):
        self.name = "precipitation"
        super().__init__()

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, fine_mask_len: int):

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)

        logMhalo = np.log10(mass)
        Rvirkpc = 1000 * comoving_radius
        XRvir = 2 # XRvir is how many virial radii to go out
        #Zmetal is the metalicity; tratcrit is the coolin crieteria -- both of these don't need to change
        Z_METAL = 0.3
        TRATCRIT = 10

        fbaryon=0.2 #put in correct value TODO
    #     Rvirkpc = 300 #put in correct value here
        
        fitarray = np.array([[350, 8e12, 10, 0.5, 2.7,  0.73, 1.2e-1, 1.2, 3.8e-4,  2.1], \
            [350, 8e12, 10, 0.3, 2.4,  0.74, 1.5e-1, 1.2, 4.2e-4, 2.1], \
        # [350, 8e12, 10 &  $Z_{\rm grad}$ & 3.6  &  0.68 &  $8.6 \times 10^{-2}$ & 1.1 & $4.1 \times 10^{-4}$ & 2.0 \\
        # [300, 5.1e12,  & 10 &  1.0 & 3.3  &  0.71 & 5.6 \times 10^{-2}$  & 1.2 & $1.7 \times 10^{-4}$ & 2.1 \\
            [300, 5.1e12,  10,  0.5, 2.6,  0.72, 8.0e-2, 1.2, 2.3e-4,  2.1], \
            [300, 5.1e12,  10,  0.3, 2.3, 0.73,  9.6e-2, 1.2, 2.6e-4,  2.1], \
        #300 &  $5.1 \times 10^{12}$  & 10 &  $Z_{\rm grad}$ & 3.4  &  0.67 &  $5.8 \times 10^{-2}$ & 1.1 & $2.6 \times 10^{-4}$ & 2.1 \\
        #[300, 5.1e12,  20, 1.0,  5.2, 0.70, 2.9e-2, 1.1, 1.0e-4, 2.1],
            [300, 5.1e12, 20, 0.5,  4.1,  0.71, 4.3e-2, 1.1, 1.4e-4, 2.1], \
            [300, 5.1e12, 20, 0.3, 3.6,  0.71, 5.1e-2, 1.2, 1.6e-4, 2.1], \
        #300 &  $5.1 \times 10^{12}$  & 20 &  $Z_{\rm grad}$ & 5.4  &  0.65 &  $3.0 \times 10^{-3}$ & 1.1 & $1.6 \times 10^{-4}$ & 2.1 \\
        # 250 &  $2.9 \times 10^{12}$  & 10 &  1.0 & 3.7  &  0.71 &  $2.8 \times 10^{-2}$ & 1.2 & $8.1 \times 10^{-5}$ & 2.2 \\
            [250, 2.9e12, 10,  0.5, 2.8, 0.72, 4.2e-2, 1.2, 1.1e-4, 2.2], \
            [250, 2.9e12, 10,  0.3, 2.4,  0.72, 5.1e-3, 1.2, 1.3e-4, 2.2], \
        #250 &  $2.9 \times 10^{12}$  & 10 &  $Z_{\rm grad}$ & 3.7 &  0.66  &  $2.9 \times 10^{-2}$ & 1.1 & $1.3 \times 10^{-4}$ & 2.1 \\
        #220 &  $2.0 \times 10^{12}$  & 10 &  1.0 & 4.0  &  0.70    &  $1.7 \times 10^{-2}$ & 1.2 & $4.2 \times 10^{-5}$ & 2.3 \\
            [220, 2.0e12,  10,  0.5, 3.0,  0.71, 2.5e-2, 1.2, 6.1e-5, 2.2], \
            [220, 2.0e12,  10, 0.3, 2.6,  0.71, 3.1e-2, 1.2, 7.2e-5, 2.2], \
        #220 &  $2.0 \times 10^{12}$  & 10 &  $Z_{\rm grad}$ & 4.1  &  0.65 &  $1.8 \times 10^{-2}$ & 1.1 & $7.4 \times 10^{-5}$ & 2.2 \\
        #220 &  $2.0 \times 10^{12}$  & 20 &  1.0 & 6.3  &  0.69  &  $8.6 \times 10^{-3}$ & 1.1 & $2.3 \times 10^{-5}$ & 2.3 \\
            [220, 2.0e12,  20,  0.5, 4.8, 0.70, 1.3e-2, 1.1, 3.4e-5, 2.2], \
            [220, 2.0e12,  20,  0.3, 4.1,  0.70, 1.6e-2, 1.1, 4.2e-5, 2.2], \
        #220 &  $2.0 \times 10^{12}$  & 20 &  $Z_{\rm grad}$ & 6.5  &  0.63 &  $9.2 \times 10^{-3}$ & 1.1 & $4.3 \times 10^{-4}$ & 2.1 \\
        # 180 &  $1.1 \times 10^{12}$  & 10 &  1.0 & 5.6  &  0.71  &  $5.4 \times 10^{-3}$ & 1.2 & $9.6 \times 10^{-6}$ & 2.2 \\
            [180, 1.1e12, 10,  0.5, 4.0,  0.71, 8.7e-3, 1.2, 1.6e-5, 2.2], \
            [180, 1.1e12, 10,  0.3, 3.4,  0.71, 1.1e-2, 1.2, 2.1e-5, 2.2], \
        #180 &  $1.1 \times 10^{12}$  & 10 &  $Z_{\rm grad}$ & 5.7  &  0.63 &  $5.9 \times 10^{-3}$ & 1.1 & $2.3 \times 10^{-5}$ & 2.2 \\
            [150, 6.3e11, 10,  0.5, 6.1,  0.68, 2.8e-3, 1.2, 7.4e-6, 2.3], \
            [150, 6.3e11, 10, 0.3, 4.9,  0.68, 3.4e-3, 1.1, 9.8e-6, 2.2], \
            [150, 6.3e11, 10, 0.1, 3.0, 0.69, 8.1e-3, 1.2, 1.9e-5, 2.2], \
            [120, 3.2e11, 10, 0.5, 6.1,  0.69, 1.4e-3, 1.2, 2.3e-6, 2.2], \
            [120, 3.2e11, 10, 0.3, 5.0,  0.70, 1.9e-3, 1.2, 2.9e-6, 2.2], \
            [120, 3.2e11, 10, 0.1, 3.1,  0.71, 3.7e-3, 1.2, 5.1e-6, 2.2], \
            [120, 3.2e11, 20, 0.5, 9.6,  0.69, 7.0e-4, 1.2, 1.2e-6, 2.2],\
            [120, 3.2e11, 20, 0.3, 7.9,  0.70, 9.5e-4, 1.2, 1.5e-6, 2.2],\
            [120, 3.2e11, 20, 0.1, 4.9,  0.71, 1.9e-3, 1.2, 2.7e-6, 2.2]])

        #chooses metalicity and cooling criterion parameters we will interpolate on
        reducedarr = fitarray[(fitarray[:, 3] == Z_METAL) & (fitarray[:, 2] ==  TRATCRIT)]
        reducedarr = reducedarr[::-1] #reverses array           

        logMhalo_zp2 = logMhalo*((1+0.2)/(1+redshift))**(3/2)
        #better interpolation
        logn1 = interp1d(np.log10(reducedarr[:, 1]), np.log10(reducedarr[:, 6]), kind='linear', fill_value='extrapolate')(logMhalo_zp2)   
        n1 = 10**logn1
        xi1 = interp1d(np.log10(reducedarr[:, 1]), reducedarr[:, 7], kind='linear', fill_value='extrapolate')(logMhalo_zp2) 
        logn2 = interp1d(np.log10(reducedarr[:, 1]), np.log10(reducedarr[:, 8]), kind='linear', fill_value='extrapolate')(logMhalo_zp2)   
        n2 = 10**logn2
        xi2 = interp1d(np.log10(reducedarr[:, 1]), reducedarr[:, 9], kind='linear', fill_value='extrapolate')(logMhalo_zp2) 
            
        #print(logMhalo, n1, xi1, n2, xi2)

        rkpc = np.logspace(0, np.log10(Rvirkpc*XRvir), 300)#number of radial bins

        #Voit 2018 fitting formulae
        rhoarr = np.array([rkpc, 1/np.sqrt(1/(n1*(rkpc)**-xi1+ 1e-20)**2 + 1/(n2*(rkpc/100)**-xi2 + 1e-20)**2)])
        
        #Integrate to see how much mass is missed by this profile  (I've checked these seems reasonable)
    #     rhointerp = interp1d(np.log(rhoarr[0]), 4.*np.pi*rhoarr[0]**3*rhoarr[1], kind='cubic', fill_value='extrapolate')
        rhointerp = interp1d(np.log(rhoarr[0]), 4.*np.pi*rhoarr[0]**3*rhoarr[1], kind='linear', fill_value='extrapolate')
    #     conv = (1.989e33/(1.67e-24))/(3.086e21)**3  #constants: use your values, should have mean molecular weight, which is 1.2 
        conv = (msun*1E3/(mu*mp))/KPCTOCM**3
        conv1 = (1.989e33/(1.67e-24))/(3.086e21)**3  #constants: use your values TODO
        assert (np.isclose(conv, conv1))
        mtotal = integrate.quad(rhointerp, 0, np.log(XRvir*Rvirkpc))[0]/conv



        
        #add in rest of mass 
        neconstant =(10**logMhalo*fbaryon-mtotal)/(4.*np.pi*(XRvir*Rvirkpc)**3)*conv
    #     print(neconstant)
        length = len(rkpc)
    #     print("shape = ", length)
        for i in range(length):
            if rhoarr[0,  i] < XRvir*Rvirkpc:
    #             rhoarr[1, i] = np.array([1/np.sqrt(1/(n1*(rkpc)**-xi1+ 1e-20)**2 + 1/(n2*(rkpc/100)**-xi2 + 1e-20)**2)])[1,i] + neconstant
                rhoarr[1, i] = rhoarr[1, i] + neconstant
        
        #print("ftotal =", mtotal/10**logMhalo/fbaryon, neconstant) #.2 is fraction of baryons
        cellsize_kpc = cellsize * 1000 # kpc
        
    #     x = np.ogrid[-10*resolution: 10*resolution]
        y,x = np.ogrid[-20*resolution: 20*resolution, -20*resolution:20*resolution] # TODO delete this line

    #     f1= lambda x, y, z: my_func(((x**2+y**2+z**2)**.5), n1,n2,xi1,xi2,neconstant,cellsize_kpc)
        f1= lambda x, y, z: precipitation_func(((x**2+y**2+z**2)**.5), n1,n2,xi1,xi2,neconstant,cellsize_kpc,Rvirkpc,XRvir,redshift)
                    
        vec_integral=np.vectorize(project_spherical_3Dto2D)
        
        mask1 = vec_integral(f1,x,y,XRvir*Rvirkpc/cellsize_kpc)
        r=(x**2+y**2)**.5 # * scale_down
        mask1=mask1.astype(float)
        mask1[r > (XRvir*Rvirkpc/cellsize_kpc)]=0         

    #     mask1[r <= (XRvir*Rvirkpc/cellsize_kpc)] =+ neconstant 
            
        return mask1

# The user can create their own function
# All length scales have to be converted into units of cellsize
def precipitation_func(r, n1,n2,xi1,xi2,neconstant,cellsize_kpc,Rvirkpc,XRvir,redshift):
        
    #final_ar = np.array([1/np.sqrt(1/(n1*((r+.5)*cellsize_kpc)**-xi1)**2 + 1/(n2*((r+.5)*cellsize_kpc/100)**-xi2)**2)])# + neconstant
    x = (np.array(r) <= XRvir*Rvirkpc)

    # TODO I got this line broken from Adnan and had to guess where the ) goes to make it valid
    final_ar =   np.array([1/np.sqrt(1/(n1*((r+.5)*cellsize_kpc/(1+redshift))**-xi1)**2 + 1/(n2*((r+.5)*cellsize_kpc/(100*(1+redshift))**-xi2)**2)) + x.astype(int)*neconstant])

    return final_ar


# This function subtracts the halos from the density field

# arguments: 
# haloArray: dataframe of halos, sorted by mass.
# bin_markers: array giving the indexes of halos at the edges of mass bins
# profile: tophat, NFW etc
# scaling_radius: scale radius for tophat halos
def subtract_halos(haloArray, bin_markers, profile: CGMProfile, scaling_radius: float, redshift: float):
    
    # TODO I think this is effectively hardcoded to the 256 Bolshoi grid size.
    df = haloArray
    no_cells = 1024
    cellsize = L/(1024) 
    chunks = len(bin_markers) - 1
    
    # array of halo masses and radii
    Mvir_avg = np.zeros(chunks)
    conv_rad = np.zeros(chunks)

    # convolution mask array
    convolution = np.zeros([chunks,no_cells,no_cells])    
    
    # fine and coarse map settings
    fine_mask_len = 10
    scale_down = 2  # making the grid coarser
    nbig = fine_mask_len*2
    nsmall = int(nbig/scale_down)
    
    # loops through the list of dataframes each ordered by ascending mass
    for j in range(0,chunks):

        if bin_markers[j] == bin_markers[j+1]:
            # Special case - mass bin is empty. Just set results to 0 for this mass bin.
            Mvir_avg[j] = 0
            conv_rad[j] = 0
            convolution[j,:,:] = 0
            continue

        Mvir_avg[j] = np.mean((df['Mvir'][bin_markers[j]:bin_markers[j+1]]))/h
        conv_rad[j] = comoving_radius_for_halo(Mvir_avg[j], redshift) # comoving radius

        fine_mask = profile.get_mask(Mvir_avg[j], conv_rad[j], redshift, 1, scaling_radius, cellsize, fine_mask_len)

        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
        
        # Area of cells needed for normalization
        totalcellArea4 = sum(sum(coarse_mask))* ((cellsize)**2)

        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])   
        
        # The coordinates are being multiplied by 4 to yield the halo coordinates on the 1024 grid
        ix = ((((np.around(4*((df[bin_markers[j]:bin_markers[j+1]]['x'].values)/(250/256))))))%(1024)).astype(int)
        iy = ((((np.around(4*((df[bin_markers[j]:bin_markers[j+1]]['y'].values)/(250/256))))))%(1024)).astype(int)
        xy=(ix,iy)

        # BUG issue: the method does not add repeated coordinates
        halo_cell_pos[xy] += 1
        
        # convolve the mask and the halo positions
        convolution[j,:,:] = (Mvir_avg[j]/(totalcellArea4)) * my_convolve(halo_cell_pos,coarse_mask)    
        
    return (convolution.sum(0))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM)




# Status: This is the latest convolution function

# The function add halos

# arguments: 
# haloArray: dataframe of halos, sorted by mass.
# resolution: int of resolution, will be multiplied by 1024
# bin_markers: array giving the indexes of halos at the edges of mass bins
# profile: tophat, NFW etc
# scaling_radius: scale radius for tophat halos
def add_halos(haloArray, resolution: int, bin_markers, profile: CGMProfile, scaling_radius: int, redshift: float):
    
    df = haloArray
    no_cells = 1024 * resolution
    cellsize = L / no_cells
    chunks = len(bin_markers) - 1

    # array of halo masses and radii
    Mvir_avg = np.zeros(chunks)
    conv_rad = np.zeros(chunks)

    # convolution mask array
    convolution = np.zeros([chunks,no_cells,no_cells])
    
    # fine and coarse map settings
    fine_mask_len = 20*resolution # TODO CGMBrush is forcing this mask size scheme. Should implementers get to choose?
    scale_down = 2  # making the grid coarser    
    nbig = fine_mask_len*2
    nsmall = int(nbig/scale_down)

    # store all profile masks
    addition_masks =np.zeros([chunks,nsmall,nsmall])
    
    # loops through the list of dataframes each ordered by ascending mass
    for j in range(0,chunks):

        if bin_markers[j] == bin_markers[j+1]:
            # Special case - mass bin is empty. Just set results to 0 for this mass bin.
            Mvir_avg[j] = 0
            conv_rad[j] = 0
            addition_masks[j,:,:] = 0
            convolution[j,:,:] = 0
            continue

        Mvir_avg[j] = np.mean((df['Mvir'][bin_markers[j]:bin_markers[j+1]])) / h
        conv_rad[j] = comoving_radius_for_halo(Mvir_avg[j], redshift) # comoving radius

        fine_mask = profile.get_mask(Mvir_avg[j], conv_rad[j], redshift, resolution, scaling_radius, cellsize, fine_mask_len)

        if profile.debug:
            with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
                print("FINE MASK %s" % j)
                print(fine_mask)

        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)

        if profile.debug:
            with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
                print("COARSE MASK %s" % j)
                print(coarse_mask)

        # Area of cells needed for normalization
        totalcellArea4 = sum(sum(coarse_mask))* ((cellsize)**2)
        
        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])    
        
        # The coordinates are being multiplied by 4 to yield the halo coordinates on the 1024 grid
        # TODO the resolution math here can be written in a more intuitive way
        ix = ((((np.around(4*resolution*((df[bin_markers[j]:bin_markers[j+1]]['x'].values)/(250/256))))))%(resolution*1024)).astype(int)
        iy = ((((np.around(4*resolution*((df[bin_markers[j]:bin_markers[j+1]]['y'].values)/(250/256))))))%(resolution*1024)).astype(int)  
        xy=(ix,iy)

        # Adnan: BUG issue: the method does not add repeated coordinates. Ian: Is that right?
        halo_cell_pos[xy] += 1

        # convolve the mask and the halo positions
        convolution[j,:,:] = (Mvir_avg[j]/(totalcellArea4)) * my_convolve(halo_cell_pos,coarse_mask)

        # store addition masks
        addition_masks[j,:,:]= (Mvir_avg[j]/(totalcellArea4))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM)*coarse_mask
        
    return (convolution.sum(0))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM), conv_rad, addition_masks, Mvir_avg


# Halos removed field

#This function combines many steps of subtraction and smoothing to yield a density field from which halos have been removed
def halos_removed_field(current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,subtraction_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc):
    
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    bin_markers= halo_array_for_convolution[1]
    
    # convolve halos
    subtraction_profile = subtract_halos(df,bin_markers,subtraction_halo_profile,scaling_radius,redshift)
    assert not np.any(np.isnan(subtraction_profile))
    subtraction_profile_smooth = gauss_sinc_smoothing(subtraction_profile,sigma_gauss,width_sinc,1)
    assert not np.any(np.isnan(subtraction_profile_smooth))
    
    # create coarse grid
    subtraction_coarse= smoothfield(subtraction_profile_smooth,1024,den_grid_size)
    assert not np.any(np.isnan(subtraction_coarse))
    
    # remove halos from the density field
    halos_removed_coarse=removeConvolvedHalos(density_field,subtraction_coarse)
    assert not np.any(np.isnan(halos_removed_coarse))
    
    return halos_removed_coarse,subtraction_coarse,subtraction_profile,subtraction_profile_smooth


# Function subtracts and adds halos
def convolution_all_steps_final(current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,halos_removed_coarse,
                       addition_halo_profile: CGMProfile,scaling_radius,resolution,sigma_gauss,width_sinc):
    
    # setup inputs for convolution
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    bin_markers= halo_array_for_convolution[1]
    
    # convolve halos for adding back
    addition_profile_initial=add_halos(df,resolution,bin_markers,addition_halo_profile,scaling_radius,redshift)
    addition_profile = addition_profile_initial[0]
    addition_profile_masks=addition_profile_initial[2]
    
    # add halos to the subtracted field
    halosremoved_fine = (np.repeat((np.repeat(halos_removed_coarse,(1024/den_grid_size)*resolution,axis=0)),(1024/den_grid_size)*resolution,axis=1))
    roll_by = int((addition_profile.shape[0]/den_grid_size)/2)
    halosremoved_fine = np.roll(halosremoved_fine, -1*roll_by, axis=0)
    halosremoved_fine = np.roll(halosremoved_fine, -1*roll_by, axis=1)
    halos_added = addition_profile +  halosremoved_fine
    
    virial_rad = addition_profile_initial[1]
    halo_masses = addition_profile_initial[3]
    
    return halos_added, halosremoved_fine,addition_profile,addition_profile_masks,halos_removed_coarse,virial_rad,halo_masses
    


# Multiple Redshift Convolution
def halo_subtraction_addition(sim_provider : SimulationProvider,den_grid_size,RS_array,min_mass,max_mass,log_bins,subtraction_halo_profile,
                             addition_profile: CGMProfile,scaling_radius,resolution):

    # Details of halo profiles

    # sigma is the variance of the gaussian. sigma**2 = 16 or 3*sigma**2 = 16
    sigma_gauss = 4

    # nsinc is the full width of the box function. So 4 means two boxes on either side of the point.
    width_sinc = 4

    halos_reAdded = np.zeros([len(RS_array),1024*resolution,1024*resolution])
    
    halos_subtraction_coarse = np.zeros([len(RS_array),den_grid_size,den_grid_size])
    halo_field = np.zeros([len(RS_array),1024*resolution,1024*resolution])
    halo_masks = np.zeros([len(RS_array),int(log_bins-1),20*resolution,20*resolution]) # This should be the same as fine_mask_len in add_halos function
    halos_removed_fields = np.zeros([len(RS_array),den_grid_size,den_grid_size]) # This should be the same as fine_mask_len in add_halos function
    
    for i in range(0, len(RS_array)):
        
        redshift = RS_array[i]
        density_field = sim_provider.get_density_field(redshift, den_grid_size)
        halos = sim_provider.get_halos(redshift)

        halos_removed = halos_removed_field(halos,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,subtraction_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc)
              
        conv_all_steps = convolution_all_steps_final(halos,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,halos_removed[0],
                           addition_profile,scaling_radius,resolution,sigma_gauss*resolution,width_sinc)
        
        
        halos_reAdded[i,:,:] = conv_all_steps[0]
        
        halos_subtraction_coarse[i,:,:] = halos_removed[1]
        halo_field[i,:,:] = conv_all_steps[2]
        halo_masks[i,:,:,:]= conv_all_steps[3]
        halos_removed_fields[i,:,:] =  halos_removed[0]
        
        
    # returns halo array and masks used to add halos back
    return halos_reAdded,halo_masks,halos_subtraction_coarse,halo_field,halos_removed_fields, conv_all_steps[5],conv_all_steps[6],halos_removed[2],halos_removed[3]



# TODO clean this up, we've basically made it do nothign over halo_subtraction_addition
def hist_profile(sim_provider: SimulationProvider, den_grid_size, RS_values, min_mass, max_mass,
                                       log_bins, subtraction_halo_profile, addition_profile: CGMProfile, scaling_radius, resolution):
    """
    This function runs the convolution code (subtraction and addition).

    Outputs: halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, halos removed field, virial radii, halo masses
    """
    
    # halo array
    t = halo_subtraction_addition(sim_provider,den_grid_size,RS_values,min_mass,max_mass,
                                       log_bins,subtraction_halo_profile,addition_profile,scaling_radius,resolution)
    # Halos-readded field
    t1 = t[0]
    
    # Halo addition masks
    t2 = t[1]
    
    # Halos subtraction coarse
    t3 = t[2]
    
    # Halo addition field
    t4 = t[3]
    
    # Halos removed field
    t5 = t[4]
    
    #if len(RS_array)==1:
    #    t6 = t1
    
    #else:
    #    t6 = stack_all_arrays(t1,RS_array)
        
    #t7 = create_histograms(t6, resolution*1024)
    #t7 = (0,0,0) 
    
    t8 = t[5]
    t9 = t[6]
    #t10 = np.zeros([len(RS_array),resolution*1024,resolution*1024])

    #for i in range(0,len(RS_array)):
    #    t10[i,:,:] = redshifted_DM(t6[i,:,:], RS_array[i])

    # Outputs: 
    # halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, halos removed field, virial radii, halo masses
    return t1,t2,t3,t4,t5,t8,t9,t[7],t[8]




####################################
# DM profile for a mass bin
####################################

# Create radial profile of a field
def radial_profile(data, center_x,center_y):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 


def make_halo_square(DM_field, ix, iy, halo_index, crop_grid):
    """Createa a square cutout of the DM field from around the ix, iy point 
    for the halo_index halo with dimensions crop_grid x crop_grid."""
    
    res=DM_field.shape[0]    
    DM_rad = np.zeros([crop_grid,crop_grid])
    i = halo_index
            
    # BUG I think this ignose halos near the edge? But then we use them in the average? Hmm...
    if ix[i] > int(crop_grid) and ix[i] < res - int(crop_grid) and iy[i] > int(crop_grid) and iy[i] < res- int(crop_grid):
        trimmed = DM_field[ix[i]-int(crop_grid/2):ix[i]+int(crop_grid/2),iy[i]-int(crop_grid/2):iy[i]+int(crop_grid/2)]
        DM_rad = trimmed

    return DM_rad


# Single function for radial profile of DM for a given DM array and grid size
def DM_vs_radius(DM_field, halo_data_frame, crop_grid_dim, bin_markers):    
    
    # TODO: Matt says we can calculate using every 1/4 the pixels or something.
    
    num_bins = len(bin_markers) -1 
    center = int(crop_grid_dim/2)
    
    DM_mass_bin= np.zeros([num_bins, crop_grid_dim, crop_grid_dim]) # TODO previously this had an extra empty bin at end I think. mistake
    rad_prof_mass_bin = []

    res=DM_field.shape[0]
    ix = ((np.around((4*res/1024)*((halo_data_frame['x'].values)/(250/256))))%(res)).astype(int)
    iy = ((np.around((4*res/1024)*((halo_data_frame['y'].values)/(250/256))))%(res)).astype(int)

    # Go through each mass bin
    for i in range(0, num_bins):
        # get a square cutout from the DM field around the center of each halo within this mass bin
        num_halos = bin_markers[i+1] - bin_markers[i]
        #halo_squares = DM_for_mass_bin(DM_field, halo_data_frame[bin_markers[i]:bin_markers[i+1]], crop_grid_dim)

        # Take all the halos in this mass bin, and get the mean DM for each pixel in the trimmed square.
        # Running sum for memory reasons
        for halo_index in range(bin_markers[i], bin_markers[i+1]):
            halo_square = make_halo_square(DM_field, ix, iy, halo_index, crop_grid_dim)
            DM_mass_bin[i,:,:] = DM_mass_bin[i,:,:] + halo_square

        DM_mass_bin[i,:,:] = DM_mass_bin[i,:,:] / num_halos # finishing computing average
        rad_prof_mass_bin.append(radial_profile(DM_mass_bin[i,:,:],center,center))
    
    arr = np.array(rad_prof_mass_bin)
    return arr, DM_mass_bin
    
# Calculate profiles of masks 
def profile_of_masks(mask_array):
    
    mask_prof_ar_shape = mask_array.shape
    prof_num=mask_prof_ar_shape[0]
    mask_center=int(mask_prof_ar_shape[1]/2)
    
    dim_rad_ar=(radial_profile(mask_array[0,:,:], mask_center,mask_center)).shape[0]

    prof_masks = np.zeros([prof_num,dim_rad_ar])
    
    for i in range(0,prof_num):
        prof_masks[i,:] = radial_profile(mask_array[i,:,:], mask_center,mask_center)
        
    return prof_masks
    
# Translate arrays by random number
#def translate_array(array,seed_x,seed_y):
def translate_array(array):
    
    # Dimension of input array
    dim1 = array.shape[0]
    translated_array= np.zeros([dim1,dim1])
    
    # Random number to translate the list by
    #random.seed(seed_x)
    trans_x = random.randint(int(dim1/4),dim1)
    #random.seed(seed_y)
    trans_y = random.randint(int(dim1/4),dim1)
    
    
    for i in range(0,dim1):
        for j in range(0,dim1):
                translated_array[i,j] = array[(i+trans_x)%dim1,(j+trans_y)%dim1]
    
    return translated_array

def empty_stack(num_boxes,dim):
    empty_array = np.zeros([num_boxes,dim,dim])
    return empty_array

def complete_stacking(stack,to_redshift):
    
    dim_stack = stack.shape[0]
    
    num_arrays_to_stack=int(numBoxes(to_redshift))
    empty_array_for_stacking = empty_stack(num_arrays_to_stack,dim_stack)


    empty_array_for_stacking[0,:,:] = stack*(1+avgZ(1))
    for i in range(1,num_arrays_to_stack):
        translated_stack = translate_array(stack)*(1+avgZ(i+1))
        empty_array_for_stacking[i,:,:] = translated_stack
    


    stacked_array = sum(empty_array_for_stacking[:,:,:]) # sums the outermost array
    
    return stacked_array

def translate_field_stack(halos_reAdded, RS_array, seed):

    halos_reAdded_translated = np.zeros(halos_reAdded.shape)
    halos_reAdded_translated[0] = halos_reAdded[0] # don't translate the 1st redshift's field, so just copy over
    
    # We want to random translations for each of the blocks. But we want to be able to make the same translations from run to run
    # so different CGM profiles can be compared. Thus, set the random seed here.
    if seed is not None:
        random.seed(seed)

    # Translate all but the first z slice
    for i in range(1, len(RS_array)):
        halos_reAdded_translated[i,:,:] = translate_array(halos_reAdded[i,:,:]) 
        #halos_reAdded_translated[i,:,:] = redshifted_DM(translate_array(halos_reAdded[i,:,:]),RS_array[i])
    
    return halos_reAdded_translated

def create_histograms(halos_reAdded_translated, resolution: int):
    """Creates histograms from a field. Resolution provided should be actual full resolution."""
    nbin=80
    hist = histArray(sum(halos_reAdded_translated),nbin, resolution, 0, 3*np.mean(sum(halos_reAdded_translated)))

    return hist








###########################################
# Utils. These are to help run the code. 
# TODO Perhaps they will go in a different part of the package.
###########################################

TEST_DIR = "data" # folder under cgmbrush/test/ that holds data files used in tests is in version control
# TODO there are tests that require the Bolshoi sims files, but those are big and not in version control.
#   Possible solution - reduce the Bolshoi files down to 1/100 the size for testing only.

def saveFig(filename_base, fig, **kwargs):
    """Saves matplotlib figures to a specified folder (defaults to a var folder outside verion control)."""
    file_path = os.path.join(VAR_DIR, filename_base)
    
    if not(os.path.exists(VAR_DIR)):
        os.makedirs(VAR_DIR)

    fig.savefig(file_path, **kwargs)

def saveResults(filename, folder = VAR_DIR, **arrays):
    """Saves numpy arrays to a specified folder (defaults to a var folder outside verion control)."""
    file_path = os.path.join(folder, filename)
    
    if not(os.path.exists(folder)):
        os.makedirs(folder)

    np.savez(file_path, **arrays) 

def saveArray(filename, *arrays, folder = VAR_DIR):
    """Saves numpy arrays to a specified folder (defaults to a var folder outside verion control)."""
    file_path = os.path.join(folder, filename)
    
    if not(os.path.exists(folder)):
        os.makedirs(folder)

    if len(arrays) == 1: # unwrap it out of tuple or it will not round-trip
        np.save(file_path, arrays[0]) 
    else:
        np.save(file_path, arrays) 
    
def loadArray(filename, folder = VAR_DIR):
    """Loads numpy arrays that were saved with saveArray."""
    file_path = os.path.join(folder, filename + ".npy")
    try:
        return np.load(file_path, allow_pickle=True)
    except FileNotFoundError:
        file_path = os.path.join(folder, filename + ".txt")
        return np.load(file_path, allow_pickle=True)


class Configuration:
    """
    Configuration for a run of cgmbrush. 
    
    Contains various utilities for running the code and analyzing the results.
    """

    # Default options
    def __init__(self, addition_profile: CGMProfile, scaling_radius, provider: SimulationProvider = None, folder=VAR_DIR, resolution=1, den_grid_size=256, RS_array=[0]):
        
        # Profile to use for adding in CGM
        self.addition_profile = addition_profile
        self.scaling_radius = scaling_radius
        self.provider = provider
        self.resolution = resolution # x1024
        self.folder = folder

        # Resolution: choose between 256 and 512 grid TODO this is Bolshoi specific
        if den_grid_size != 256 and den_grid_size != 512:
            raise ValueError("Only resolutions 256 and 512 are allowed")
        self.den_grid_size = den_grid_size 

        # User provides a redshift array
        self.RS_array = RS_array # For a single box, we only use the redshift 0 box

        # Profile used for subtracting halos from the density field
        self.subtraction_halo_profile = NFWProfile()


        self.min_mass = 10**10 # halos smaller than this shouldn't have much of a CGM
        self.max_mass = 9*10**15 # this is a little bigger than the biggest halo in Bolshoi
        self.log_bins = 30
        self.datestamp = str(datetime.date.today())
        self.seed = None

        self.npz = None
        self.results = None
        self.figure = None
        self.final_field = None
        self.add_masks = None
        self.addition_field = None
        self.sub_fine_unsmoothed = None
        self.sub_fine_smoothed = None
        self.sub_coarse = None
        self.removed_field = None
        self.virial_radii = None
        self.halo_masses = None
        self.DM_vs_R1 = None
        self.mask_profiles = None

        self.stacked_npz = None
        self.stacked_orig_field = None
        self.stacked_removed_field = None
        self.stacked_addition_field = None
        self.stacked_final_field = None




    def __del__(self):
        if (self.npz is not None):
            self.npz.close()
        if (self.stacked_npz is not None):
            self.stacked_npz.close()

    def get_filename(self):
        scaling = ''
        if self.scaling_radius > 1:
            scaling = '_' + str(self.scaling_radius)
        z_str = ''
        if self.RS_array != [0]:
            z_str = '_z'
            for z in self.RS_array: # z_0.1_0.2_0.3 for example
                z_name = self.provider.get_z_name(z)
                z_str += '_' + z_name
        return self.addition_profile.name + str(self.resolution) + scaling + '_' + str(self.den_grid_size) + z_str + "_" + self.datestamp

    def convert_and_save(self):
        """Makes results available as a dictionary, which is how the .npz archive is formatted. 

        hist_profile returns a tuple. Reading .npy files gets a tuple or array. We want everything converted to dictionary and npz.
        
        Note that .npz archives are loaded lazilly, which is important for managing memory usage."""
        # hist_profile returns a tuple. Reading .npy files gets a tuple or array. 
        if isinstance(self.results, np.ndarray):
            self.results = tuple(self.results)
        if type(self.results) is tuple:
            self.results = { 'final_density_field': self.results[0], 'add_masks': self.results[1], 'sub_coarse': self.results[2], 'add_density_field': self.results[3], 'removed_density_field': self.results[4], 'vir_radii': self.results[5], 'halo_masses': self.results[6], 'sub_fine_unsmoothed': self.results[7], 'sub_fine_smoothed': self.results[8] }
            saveResults(self.get_filename(), **self.results, folder=self.folder)
            self.final_field = self.results['final_density_field']
            self.add_masks = self.results['add_masks']
            self.sub_coarse = self.results['sub_coarse']
            self.addition_field = self.results['add_density_field']
            self.removed_field = self.results['removed_density_field']
            self.virial_radii = self.results['vir_radii']
            self.halo_masses = self.results['halo_masses']
            self.sub_fine_unsmoothed = self.results['sub_fine_unsmoothed']
            self.sub_fine_smoothed = self.results['sub_fine_smoothed']

        if type(self.results) is not dict:
            raise ValueError('Results are in an unexpected format: %s' % type(self.results))

    def get_final_field(self):
        if self.final_field is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.final_field = self.npz.get('final_density_field')
        
        return self.final_field

    def get_addition_masks(self):
        if self.add_masks is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.add_masks = self.npz.get('add_masks')
        
        return self.add_masks

    def get_subtraction_coarse_field(self):
        if self.sub_coarse is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.sub_coarse = self.npz.get('sub_coarse')
        
        return self.sub_coarse

    def get_addition_field(self):
        if self.addition_field is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.addition_field = self.npz.get('add_density_field')
        
        return self.addition_field

    def get_removed_field(self):
        if self.removed_field is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.removed_field = self.npz.get('removed_density_field')
        
        return self.removed_field

    def get_virial_radii(self):
        if self.virial_radii is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.virial_radii = self.npz.get('vir_radii')
        
        return self.virial_radii

    def get_halo_masses(self):
        if self.halo_masses is None:
            
            if self.npz is None:
                self.run(load_from_files = True)
            
            self.halo_masses = self.npz.get('halo_masses')
        
        return self.halo_masses

    def run(self, plots=False, trace=False, results_in_memory=True, load_from_files=False):
        """Run this configuration."""

        filename = self.get_filename()
        file_path = os.path.join(self.folder, filename + ".npz")

        if load_from_files:
            try:
                self.npz = np.load(file_path, allow_pickle=True)

                # This is the non-lazy approach
                #for file in npz.files:
                #    self.results[file] = npz[file] 
                #npz.close()
                            
            except IOError:
                print("Cache miss: " + filename)
                #pass # File cache doesn't exist, swallow and compute it instead

        if self.npz is None:
            print("Performing Calculations for " + filename)
                                   
            if trace:
                pr = cProfile.Profile()
                pr.enable()

            self.results = hist_profile(self.provider, self.den_grid_size, self.RS_array, self.min_mass, 
                                                self.max_mass, self.log_bins, self.subtraction_halo_profile, 
                                                self.addition_profile, self.scaling_radius, self.resolution)
            
            self.convert_and_save()
            self.npz = np.load(file_path, allow_pickle=True)

            if trace:
                pr.disable()
                s = io.StringIO()
                ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
                ps.print_stats()

                version = 1
                # TODO auto version incrementing
                perf_file_path = os.path.join(VAR_DIR, 'perf_convo_' + filename + '_v%s.txt' % version)
                with open(perf_file_path, 'w') as f:
                    f.write(s.getvalue())

        if plots: # TODO delete?
            original = self.provider.get_density_field(0, self.den_grid_size)
            background_dm = self.get_removed_field()[0]
            cgm_only = self.get_addition_field()[0]
            density_final = self.get_final_field()[0]

            vmin = 10 # for a log color plot
            vmax = max(np.max(original), np.max(background_dm), np.max(cgm_only), np.max(density_final))
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)

            fig, axes = plt.subplots(2,2,figsize=(24, 24))
            pos = axes[0][0].imshow(original, norm=norm) 
            #fig.colorbar(pos, ax=axes[0][0])
            axes[0][0].title.set_text('Original Density Field')
            pos = axes[0][1].imshow(background_dm, norm=norm) 
            #fig.colorbar(pos, ax=axes[0][1])
            axes[0][1].title.set_text('Density minus halos')
            pos = axes[1][0].imshow(cgm_only, norm=norm) 
            #fig.colorbar(pos, ax=axes[1][0])
            axes[1][0].title.set_text('CGM Profile to add')
            pos = axes[1][1].imshow(density_final, norm=norm) 
            axes[1][1].title.set_text('Final Density Field')

            fig.colorbar(pos, ax=axes, shrink=0.85)


            self.figure = fig
            saveFig(filename + '_images', fig, folder=self.folder)
        
        if not results_in_memory:
            self.clear_results()
    

    def get_stacked_orig_field(self):

        if self.stacked_orig_field is None:

            if self.stacked_npz is None:
                self.generate_stacked_fields(load_from_files=True)
            
            self.stacked_orig_field = self.stacked_npz.get('stacked_orig_field')

        return self.stacked_orig_field

    def get_stacked_removed_field(self):

        if self.stacked_removed_field is None:

            if self.stacked_npz is None:
                self.generate_stacked_fields(load_from_files=True)
            
            self.stacked_removed_field = self.stacked_npz.get('stacked_removed_field')

        return self.stacked_removed_field

    def get_stacked_addition_field(self):

        if self.stacked_addition_field is None:

            if self.stacked_npz is None:
                self.generate_stacked_fields(load_from_files=True)
            
            self.stacked_addition_field = self.stacked_npz.get('stacked_addition_field')

        return self.stacked_addition_field

    def get_stacked_final_field(self):

        if self.stacked_final_field is None:

            if self.stacked_npz is None:
                self.generate_stacked_fields(load_from_files=True)
            
            self.stacked_final_field = self.stacked_npz.get('stacked_final_field')

        return self.stacked_final_field


    def generate_stacked_fields(self, results_in_memory=True, load_from_files=False):

        # TODO do this operation for the original, removed, and addition fields too


        #translated_file = self.get_filename() + "_translated"
        stacked_file = self.get_filename() + "_stacked"
        file_path = os.path.join(self.folder, stacked_file + ".npz")
    
        if load_from_files:
            try:
                print("Loading stacked fields... ", end="")
                self.stacked_npz = np.load(file_path, allow_pickle=True)
                print("done")        
            except IOError:
                print("Cache miss: " + stacked_file)

        if self.stacked_npz is None:

            if len(self.RS_array) > 1:
                print("Creating Stacked Fields... ", end="")

                all_orig_fields = np.zeros((len(RS_array_gen), self.den_grid_size, self.den_grid_size))
                for i in range(0, len(self.RS_array)):
                    all_orig_fields[i] = self.provider.get_density_field(self.RS_array[i])
                translated_field = translate_field_stack(all_orig_fields, self.RS_array, self.seed)
                self.stacked_removed_field = np.zeros(translated_field.shape)
                for i in range(0, len(self.RS_array)):
                    self.stacked_removed_field[i,:,:] = redshifted_DM(translated_field[i,:,:], self.RS_array[i])

                translated_field = translate_field_stack(self.get_removed_field(), self.RS_array, self.seed)
                self.stacked_removed_field = np.zeros(translated_field.shape)
                for i in range(0, len(self.RS_array)):
                    self.stacked_removed_field[i,:,:] = redshifted_DM(translated_field[i,:,:], self.RS_array[i])

                translated_field = translate_field_stack(self.get_addition_field(), self.RS_array, self.seed)
                self.stacked_addition_field = np.zeros(translated_field.shape)
                for i in range(0, len(self.RS_array)):
                    self.stacked_addition_field[i,:,:] = redshifted_DM(translated_field[i,:,:], self.RS_array[i])
                                
                translated_field = translate_field_stack(self.get_final_field(), self.RS_array, self.seed)
                self.stacked_final_field = np.zeros(translated_field.shape)
                for i in range(0, len(self.RS_array)):
                    self.stacked_final_field[i,:,:] = redshifted_DM(translated_field[i,:,:], self.RS_array[i])
    
                stacked_fields = { 'stacked_orig_field': self.stacked_orig_field, 'stacked_removed_field': self.stacked_removed_field, 'stacked_addition_field': self.stacked_addition_field, 'stacked_final_field': self.stacked_final_field }
                saveResults(stacked_file, **stacked_fields, folder=self.folder)
                
                print("done")
            else:
                raise ValueError('Generating a stacked field is only applicable with data from multiple redshifts.')

        if not results_in_memory:
            self.translated_field = None
            self.stacked_field = None
        

    def generate_DM_vs_radius_profile(self, load_from_files=False):
        profile_file = '%s_DMvsR_prof' % self.get_filename()

        if load_from_files:
            try:
                self.DM_vs_R1 = loadArray(profile_file, folder=self.folder)            
            except IOError:
                print("Cache miss: " + profile_file)

        if self.DM_vs_R1 is None:       
            print("Generating DM vs R profile")
            df = create_halo_array_for_convolution(self.provider.get_halos(self.RS_array[0]), self.min_mass, self.max_mass, self.log_bins)

            trim_dim = int(10*self.resolution)

            self.DM_vs_R1 = DM_vs_radius(self.get_final_field()[0,:,:], df[0], trim_dim, df[1]) [0]
            saveArray(profile_file, self.DM_vs_R1, folder=self.folder)
    
    def generate_profile_of_masks(self, load_from_files=False):
        mask_file = '%s_masks' % self.get_filename()

        if load_from_files:
            try:
                self.mask_profiles = loadArray(mask_file, folder=self.folder)            
            except IOError:
                print("Cache miss: " + mask_file)

        if self.mask_profiles is None:
            # using 0'th index (first available redshift) here
            print("Generating Mask Profiles")

            # zoom in on middle of masks
            full_mask_len = self.resolution * 20
            zoom_start = full_mask_len // 4
            zoom_end = 3 * zoom_start
            self.mask_profiles = profile_of_masks(self.get_addition_masks()[0, :, zoom_start:zoom_end, zoom_start:zoom_end])
            saveArray('%s_masks' % self.get_filename(), self.mask_profiles, folder=self.folder)

    def clear_results(self):
        """Clears memory-intense results from memory. Saved files are preserved. Results can be recovered quickly by running with load_from_files=True."""
        self.results = None
        self.final_field = None
        self.addition_field = None
        self.translated_field = None
        self.stacked_field = None
        # TODO is calling del(...) better? 
        #gc.collect()
        # should allow garbage collection to happen
        
