###################################################################################################
#
# cgmbrush.py 	        (c) Ian Williams, Adnan Khan, Carter Archuleta, Owen Fairbairn, Matt McQuinn
#     				    	ianw89@live.com, ianw89@uw.edu
#
###################################################################################################

from __future__ import print_function 
from __future__ import division
from distutils import core
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
from cgmbrush.cosmology import cosmology as cosmo
from cgmbrush.cosmology import halo as halo

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
# Miscelaneous functions
########################################

# Dispersion measure: analytical function
def DM_analytical(z):
    return lightspeed* integrate.romberg(integrandDM,0,z) #integrates a function with given variables

# integrand of dispersion measure
def integrandDM(z):
    return PcinMpc*cosmo.dConfDistdz(z)*cosmo.elecD(z)/(1+z)**2 

#DM is redshifted by a factor of 1+z
def redshifted_DM(DM,z):
    return DM *(1+z)  # Psc cm-2

# finding redshift for comoving distance
def RSvsCDar(z):
    RSvsCD = np.zeros((100,2))
    ar=np.linspace(0,z,100)
    for i in range(0,len(ar)):
        RSvsCD[i,0] = ar[i]
        RSvsCD[i,1] = cosmo.Dconf(ar[i])
    return RSvsCD

# function to get redshift from the comoving distance by interpolating upto z, returns Mpc/h 
def CDtoRS(CD,z):
    x = RSvsCDar(z)[:,0]
    y = RSvsCDar(z)[:,1]
    p = P.fit(x, y, 2)
    return ((p - CD).roots()[0])
    
# number of boxes needed for a given redshift
def numBoxes(z, Lbox):
    return float(round(cosmo.Dconf(z)/(Lbox)))
    
def z_eff(zmin,zmax,Lbox):
    return ((DM_analytical(zmax)-DM_analytical(zmin))/((Lbox*PcinMpc)*cosmo.elecD(0)))**(1) - 1

def RS_array_gen(z_max,Lbox):
    """Computes redshift values for each box when stacking boxes of length L out to redshift z_max.
    
    Returns an array of these z values.
    """
    RS_array = []
    RS_array.append(CDtoRS(Lbox,1))
    for i in range(1,int(numBoxes(z_max, Lbox))):
        RS_array.append(z_eff(CDtoRS(float(Lbox*i),z_max),CDtoRS(float(Lbox*(i+1)),z_max),Lbox))

    return RS_array

# normalize DM  #Matt -- Understand this beetter
def normDM(DM, z, resolution, Lbox):
  return DM *PcinMpc* Lbox/resolution * cosmo.elecD(z) /(1+z)**2  # Psc cm-2  


def DM_statistics(field):
    """
    Prints off statistics on the DM distribution of the provided DM field for analysis purposes.
    """
    per = [50,80,95,99,99.9,99.99,99.999,99.9999]

    with np.printoptions(precision=2):

        print("DM FIELD STATISTICS")
        print("Number of pixels: ", np.size(field))
        print("Min: ", np.min(field))
        print("Max: ", np.max(field))

        #print("200 bin histogram")
        #hist, bins = np.histogram(field, bins=200)
        #print(hist)
        #print("Bins: ", bins)

        print("Percentiles: ", per)
        print(np.percentile(field, per))

        print("Standard Deviations")
        data = field
        std = np.round(data.std(),2)
        print('Raw:        ', std)
        data = data[0:len(data)-1]
        print('Drop 1:     ', np.round(data.std(),2))
        data = data[0:len(data)-9]
        print('Drop 10:    ', np.round(data.std(),2))
        data = data[0:len(data)-90]
        print('Drop 100:   ', np.round(data.std(),2))
        #data = data[0:len(data)-900]
        #print('Drop 1000:  ', np.round(data.std(),2))
        #data = data[0:len(data)-9000]
        #print('Drop 10000: ', np.round(data.std(),2))

        print("\n\n- - - - - - - - - - - - - - - - - - - -\n")

        return std


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

    @property
    @abc.abstractmethod
    def Lbox(self):
        """Length of a boxside in Mpc/h."""
        pass

    @property
    @abc.abstractmethod
    def halofieldresolution(self):
        """Resolution of the halo field (in pixels), which may differ from the raw density field."""
        pass
    
    @abc.abstractmethod
    def get_density_field(self, redshift: float, resolution: int):
        """Gets the 2D density field (in pixels) for the given redshift and grid resolution."""
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

    Lbox = 250/cosmo.h # comoving length of Bolshoi box in Mpc
    halofieldresolution = 1024

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

        #print("Density Field Specs: ", redshift, resolution, self.Lbox, cosmo.elecD(redshift))
        
        
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
            print("Reading from raw Bolshoi Files")
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
            # x,y,z are in Mpc/h
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
        """Imports the density field for a given redshift and smooths it. Stores the result in a numpy file that is faster to access."""
        
        if resolution != 256 and resolution != 512:
            raise ValueError("Only resolution 256 or 512 is supported.")

        resStr = str(resolution)

  
        # Have to hardcode z=0 table because of the unique column names
        if z_name == '0':
            # reading density field and halos data
            file_path= os.path.join(SIMS_DIR, 'dens'+resStr+'-z-0.0.csv.gz')
            pdDens=pd.read_csv(file_path)

            # extracting columns
            pdDensN=pdDens[['Bolshoi__Dens'+resStr+'_z0__ix','Bolshoi__Dens'+resStr+'_z0__iy','Bolshoi__Dens'+resStr+'_z0__iz','Bolshoi__Dens'+resStr+'_z0__dens']]

            # 3D density array
            pdDensN=pdDensN.sort_values(['Bolshoi__Dens'+resStr+'_z0__ix','Bolshoi__Dens'+resStr+'_z0__iy','Bolshoi__Dens'+resStr+'_z0__iz'])
            tden = pdDensN['Bolshoi__Dens'+resStr+'_z0__dens'].values
            tden2 = np.reshape(tden,(resolution,resolution,resolution))
            return normDM((tden2+1).sum(2), 0, resolution, self.Lbox)

        else:
            file_path = os.path.join(SIMS_DIR, 'dens'+resStr+'-z-{}.csv.gz'.format(z_name))
            den = pd.read_csv(file_path)
            den2=den[['Bolshoi__Dens'+resStr+'__ix','Bolshoi__Dens'+resStr+'__iy','Bolshoi__Dens'+resStr+'__iz','Bolshoi__Dens'+resStr+'__dens']]

            # 3D density array
            den_sorted=den2.sort_values(['Bolshoi__Dens'+resStr+'__ix','Bolshoi__Dens'+resStr+'__iy','Bolshoi__Dens'+resStr+'__iz'])
            den_vals = den_sorted['Bolshoi__Dens'+resStr+'__dens'].values
            den = np.reshape(den_vals,(resolution,resolution,resolution))
            
            den = normDM((den+1).sum(2),0, resolution, self.Lbox) 

            test_l= (np.repeat((np.repeat(den,4,axis=0)),4,axis=1))

            test_sm =  gauss_sinc_smoothing(test_l,4,4,1, self.halofieldresolution)

            assert resolution == 256
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
#returns histogram and bins
def  create_halo_array_for_convolution(pdHalos, M_min, M_max, logchunks):

    """Creates an array of indexes corresponding to indexes in pdHalos at the edges of logarithmic mass bins.

    The length will be logchunks, and the effective number of mass bins is [logchunks-1]"""
    halos = haloArray_minmax(pdHalos,M_min,M_max)
    df = halos.sort_values(by='Mvir',ascending=True)
    sorted_haloMasses = df['Mvir'].values
    histH, binsH = np.histogram(sorted_haloMasses,bins=np.logspace(np.log10(M_min),np.log10(M_max), logchunks))
    bins=np.append([0],histH.cumsum()).tolist()


    #   #gives the mean mass in each bin
    #    mean_mass_bins = np.zeros(len(binsH)-1)
    #    for i in range(len(binsH)-1):
    #        lbool = (sorted_haloMasses[:] > binsH[i]) & (sorted_haloMasses[:] < binsH[i+1])
    #       #print("shape", np.shape(lbool),binsH[i], binsH[i+1] )
    #        mean_mass_bins[i] = np.mean(sorted_haloMasses[lbool])

    
    return df, bins

def project_spherical_3Dto2D(f, x, y, cutoff, cache):
    #print("Rvir, x, y: {}, {}, {}".format(Rvir, x, y))
    boundary = math.sqrt(max(0.0, cutoff**2-(x**2+y**2)))
    if boundary == 0.0:
        return 0.0
    else:
        polar_radius = round(x**2 + y**2, 6)
        if (polar_radius in cache):
            return cache[polar_radius]
        else:
            #print("Boundary: {}".format(boundary))
            result = 2 * integrate.quad(f, 0, boundary, args=(x,y))[0] 
            cache[polar_radius] = result
            return result

def project_spherical_3Dto2D_optimized(f, x, y, cutoff):
    """Project a spherically symmetric profile from 3D to 2D.
    
    cutoff should be in whatever units x and y are in, and is the limit of integration (a hard boundary of the sphere)""" 
    vectorized_func = np.vectorize(project_spherical_3Dto2D)
    cache = {}
    #print("Projecting...")
    result = vectorized_func(f, x, y, cutoff, cache)
    #print(len(cache))
    #print(cache.keys())
    return result


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
#eta is the enhancement in resolution over the halo field
def gauss_sinc_smoothing(smoothing_array,sigma_gauss,width_sinc,eta, halofieldresolution):
    
    func = smoothing_array
    fft_func = np.fft.fft2(func)

    # sigma is the variance of the gaussian. sigma**2 = 16 or 3*sigma**2 = 16
    sigma = sigma_gauss

    # nsinc is the full width of the box function. So 4 means two boxes on either side of the point.
    nsinc = width_sinc

    # wavenumbers for the transforms
    kx = (np.fft.fftfreq(halofieldresolution*eta))
    ky = kx[:,np.newaxis]


    prod = np.zeros([halofieldresolution*eta,halofieldresolution*eta],dtype=complex)

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
        self.name = ''
        self.pretty_name = ''
        self.parameter_name_str = ''

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_mask') and
        callable(subclass.get_mask) or 
        NotImplemented)

    @abc.abstractmethod
    def get_mask(self, mass: float, comoving_rvir: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):
        """Constructs and returns a 2D mask to convolve with the halo locations
           for the specified halo parameters and redshift."""
        raise NotImplementedError

class MassDependentProfile(CGMProfile):
    """Profile to use when you want to use different profiles for different mass ranges."""

    def __init__(self, lower_profile: CGMProfile, higher_profile: CGMProfile, cutoff_mass: float):
        super().__init__()
        self.name = lower_profile.name + "_and_" + higher_profile.name + "_{:.1f}".format(math.log10(cutoff_mass))
        self.lower_profile = lower_profile
        self.higher_profile = higher_profile
        self.cutoff_mass = cutoff_mass

    def get_mask(self, mass: float, comoving_rvir: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):
        if mass <= self.cutoff_mass:
            return self.lower_profile.get_mask(mass, comoving_rvir, redshift, resolution, cellsize, fine_mask_len)
        else:
            return self.higher_profile.get_mask(mass, comoving_rvir, redshift, resolution, cellsize, fine_mask_len)


class TophatProfile(CGMProfile):

    def __init__(self, rvir_factor=1):
        super().__init__()
        self.name = "tophat"
        self.pretty_name = "Tophat"
        self.rvir_factor = rvir_factor # how many rvir's should the tophat extend to
        if (rvir_factor > 1):
            self.parameter_name_str = str(rvir_factor)

    def get_mask(self, mass: float, comoving_rvir: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):

        # fine mask
        # fine mask size has to correspond to the size of the mask that I eventually trim
        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)
        
        r = (x**2+y**2)**.5
        fine_mask = r <= (self.rvir_factor * comoving_rvir / cellsize) 
        return fine_mask.astype(float)

class SphericalTophatProfile(CGMProfile):

    def __init__(self, rvir_factor=1):
        super().__init__()
        self.name = "STH"
        self.pretty_name = "3D Tophat"
        self.rvir_factor = rvir_factor # how many rvir's should the tophat extend to
        if (rvir_factor > 1):
            self.parameter_name_str = str(rvir_factor)

    def get_mask(self, mass: float, comoving_rvir: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)

        r = (x**2+y**2)**.5

        fine_mask = r <= (self.rvir_factor * comoving_rvir / cellsize)
        fine_mask=fine_mask.astype(float)
        
        Rv = (self.rvir_factor * comoving_rvir / cellsize)
        fine_mask = fine_mask * ((1-((r/Rv)**2))**(2))**(1./4.)
        return fine_mask

class NFWProfile(CGMProfile):
    
    def __init__(self):
        super().__init__()
        #print("Initialized NFW Profile")
        self.name = "NFW"
        self.pretty_name = "NFW"
        self.vec_NFW2D = np.vectorize(self.NFW2D)

    def get_analytic_profile(self, mass: float, redshift: float):
        
        comoving_rvir = halo.comoving_rvir(cosmo, mass, redshift) # comoving radius

        
        rvals = np.logspace(np.log10(.01*comoving_rvir), np.log10(comoving_rvir, num=500, base=10))
        R_s= comoving_rvir/(halo.halo_conc(cosmo, redshift,mass))  
        rho_nought = halo.rho_0(cosmo, redshift,mass,R_s)
         
        return rvals, halo.NFW_profile(rvals,rho_nought,R_s, 0)

        
    def NFW2D(self, x, y, rho_nought, R_s, Rvir):
        """Project NFW profile from 3D to 2D.
        Rvir is in units of virial radii per cell.""" 
        offset=float(.5) 
        boundary = math.sqrt(max(0.0, Rvir**2-(x**2+y**2)))
        if boundary == 0.0: #x, y lie outside of virial radius
            return 0.0
        else:
            #return 2 * integrate.quad(lambda x, y, z: rho_nought/(((offset+(x**2+y**2+z**2)**.5)/R_s)*(1+((x**2+y**2+z**2)**.5)/R_s)**2), 0, boundary, args=(x,y))[0]
            return 2 * integrate.quad(lambda x, y, z:  halo.NFW_profile(np.sqrt(x**2+y**2+z**2),rho_nought,R_s, offset), 0, boundary, args=(x,y))[0]



    def get_mask(self, mass: float, comoving_rvir: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):
        #print("Making NFW mask for M = ", mass)
        # TODO ? add redshift to the function above ?
        R_s= comoving_rvir/(halo.halo_conc(cosmo, redshift,mass)*cellsize)  
        rho_nought = halo.rho_0(cosmo, redshift,mass,R_s)

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)
        
        if self.debug:
            with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
                print('R_s: {}, rho_0: {}, r/size: {}'.format(R_s, rho_nought, comoving_rvir / cellsize))

        # Integration bound is at the Rvir exactly because this is how the mass is defined
        integration_bound = comoving_rvir / cellsize
        fine_mask = self.vec_NFW2D(x, y, rho_nought, R_s, integration_bound)

        r=np.sqrt(x*x+y*y) 
        
        fine_mask = fine_mask.astype(float)
        fine_mask[r > comoving_rvir / cellsize] = 0 # sets very small numbers to 0

        return fine_mask

class FireProfile(CGMProfile):
    """
    This is FIRE simulation profile (from https://arxiv.org/pdf/1811.11753.pdf and using their $r^{-2}$ scaling)
    I'm using the $r^{-2}$ profile they find, with an exponetial cutoff at rmax, where we use conservation of mass to determine rmax. 
    
    So: $$\\rho = \\rho_{0}\left(\\frac{r}{r_*} \\right)^{-2} \exp[-r/r_{\\rm max}]$$

    Their results technically hold for $10^{10} - 10^{12} M_\odot$ and $z=0$, but I think it's reasonable to assume they extrapolate to other moderate redshifts and somewhat more massive halos.

    To normalize things we are using that:
    $$M_{\rm gas} = \int 4 \pi r^2 dr \rho_{0} \left(\frac{r}{r_*} \right)^{-2} \exp[-r/r_{\rm max}] = 4 \pi r_{\rm max} \rho_0 r_*^2$$

    Note that sometimes I'm working with number densities and sometimes mass densities.
    """
    
    def __init__(self):
        super().__init__()
        #print("Initialized Fire Profile")
        self.name = "fire"
        self.pretty_name = "FIRE"

    def get_mask(self, mass: float, comoving_rvir: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)             
         
        rmax,Rinterp,rho0 = self.get_Fire_params(mass, comoving_rvir, redshift) 
        #with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
        #    print("Mass {:.1e}".format(mass))
        #    print("Fire parameters: rmax={} Mpc or {} cells, Rinterp={} Mpc, rho0={}".format(rmax,rmax/cellsize,Rinterp,rho0))

        ## creating a mask        
        # x, y, z are cells, the computed r is therefor also in cells. 
        # rmax and Rinterp are in Mpc as per above, so divide to convert to cells. Epsilon is half a cellsize.
        epsilon_cells = 0.5
        fire_integral = lambda x, y, z: self.fire_func(((x**2+y**2+z**2)**.5), rmax/cellsize, Rinterp/cellsize, rho0, epsilon_cells)

        integration_bound = 4*rmax/cellsize # The exponential cutoff introduces a factor of ~0.01 by this point, hard cutoff here now.
        fine_mask = project_spherical_3Dto2D_optimized(fire_integral, x, y, integration_bound) 
        fine_mask = fine_mask.astype(float)

        return fine_mask

    # For testing:  outputs array showing analytic function for fire profile 
    def get_analytic_profile(self, mass: float, redshift: float):

        comoving_rvir = halo.comoving_rvir(cosmo, mass, redshift) # comoving radius

        rmax,Rinterp,rho0 = self.get_Fire_params(mass, comoving_rvir, redshift)

        #print("Fire parameters", rmax,Rinterp,rho0 )
        
        rvals = np.logspace(np.log10(.01*comoving_rvir), np.log10(5*rmax), num=500, base=10)
        radial_profile = self.fire_func(rvals, rmax, Rinterp, rho0, 0)
        return rvals, radial_profile

    def get_Fire_params(self, mass: float, comoving_rvir: float, redshift: float):
    
        # These are specifications taken from the fire simulations
        RinterpinRvir = 0.3  # this is the point where I read off the density nrmalization

        #profile taken from https://arxiv.org/pdf/1811.11753.pdf as described in CGM brush paper
        MfbMh = np.array([0.1,0.2,0.3,0.5,0.8,1])
        Mh = msun*np.array([10**10,10**11,10**12,10**13,10**14,10**15]) # halo masses in grams
        rv = Mpc*np.array([halo.comoving_rvir(cosmo, 10**10,0), halo.comoving_rvir(cosmo, 10**11,0), halo.comoving_rvir(cosmo, 10**12,0), halo.comoving_rvir(cosmo, 10**13,0), halo.comoving_rvir(cosmo, 10**14,0), halo.comoving_rvir(cosmo, 10**15,0)]) # radii for the given mass bins
        
        r0= .3*rv
        nHarr = MfbMh*Mh*cosmo.fb /((np.pi*4*r0**2*rv*mean_molecular_weight_electrons*mprot)) # g / cm^3
        #print(nHarr) # We give these in a footnote in the FIRE description in CGM brush paper

        nHinterp = interp1d(np.array([10., 11., 12., 13., 14., 15.]), nHarr, fill_value="extrapolate")
        if self.debug:
            print("Mass: {}, Rvir: {} Mpc, Gas Density {} g/cm^3".format(np.log10(mass), comoving_rvir, nHinterp(np.log10(mass))))

        rho0 = nHinterp(np.log10(mass)) #pivot density at r0
        Rinterp = RinterpinRvir * comoving_rvir #rvir(mass, z)
        
        Ntot = mass*msun*cosmo.fb/(mean_molecular_weight_electrons*mprot)/(Mpc**3) # number density (cm^-3)
        rmax = Ntot/(4.*np.pi*rho0*Rinterp**2 )  #from integrating above expression for rho

        # If a 3D function is desired, these two lines can be used to plot it
        #rarr = np.logspace(-2, np.log10(5*rmax), 100) # Cut off at 5 times exponential cutoff
        #rhoarr = rho0*(rarr/Rinterp)**-2*np.exp(-rarr/rmax) #number density: per cm3
        
        return rmax, Rinterp, rho0


    # epsilon is to avoid divergence from center cell
    def fire_func(self, r, rmax, Rinterp, rho0, epsilon):
        """r is the radius we are computing at (in cells)
        rmax is a parameter (in cells)
        Rinterp is another parameter (in cells)
        rho0 is a density baseline; this won't matter later because we normalize in add_halos()
        epsilon is an offset (in cells)"""
        result =  rho0 * np.exp(-r / rmax) * ((r + epsilon) / Rinterp)**-2  #Matt: the point 5 is to eliminate divergence in center, but before this was evalulated in cell sizes
        #print("Fire Projection: r={} cells rmax={} cells, Rinterp={} cells, rho0={}, result={}".format(r,rmax,Rinterp,rho0,result))
        return result


#This is the precipitation model of Voit et al (2018); https://arxiv.org/pdf/1811.04976.pdf
#as well as many other papers
class PrecipitationProfile(CGMProfile):
    """
    Voit Perciptation limited model from Appendix A in https://arxiv.org/pdf/1811.04976.pdf.
    """
    
    def __init__(self, XRvir=3, Z_METAL=0.3, epsilon=0.5):
        super().__init__()
        self.name = "precipitation"
        self.pretty_name = "Precipitation"
        self.epsilon = epsilon

        #Table taken from appendix in Voit et al (2018); https://arxiv.org/pdf/1811.04976.pdf; different entries vary metalicity; for z=0 but the above mass mapping corrects for this
        self.fitarray = np.array([
            [350, 8e12, 10, 0.5, 2.7,  0.73, 1.2e-1, 1.2, 3.8e-4,  2.1], \
            [350, 8e12, 10, 0.3, 2.4,  0.74, 1.5e-1, 1.2, 4.2e-4, 2.1], \
            [300, 5.1e12,  10,  0.5, 2.6,  0.72, 8.0e-2, 1.2, 2.3e-4,  2.1], \
            [300, 5.1e12,  10,  0.3, 2.3, 0.73,  9.6e-2, 1.2, 2.6e-4,  2.1], \
            [300, 5.1e12, 20, 0.5,  4.1,  0.71, 4.3e-2, 1.1, 1.4e-4, 2.1], \
            [300, 5.1e12, 20, 0.3, 3.6,  0.71, 5.1e-2, 1.2, 1.6e-4, 2.1], \
            [250, 2.9e12, 10,  0.5, 2.8, 0.72, 4.2e-2, 1.2, 1.1e-4, 2.2], \
            [250, 2.9e12, 10,  0.3, 2.4,  0.72, 5.1e-3, 1.2, 1.3e-4, 2.2], \
            [220, 2.0e12,  10,  0.5, 3.0,  0.71, 2.5e-2, 1.2, 6.1e-5, 2.2], \
            [220, 2.0e12,  10, 0.3, 2.6,  0.71, 3.1e-2, 1.2, 7.2e-5, 2.2], \
            [220, 2.0e12,  20,  0.5, 4.8, 0.70, 1.3e-2, 1.1, 3.4e-5, 2.2], \
            [220, 2.0e12,  20,  0.3, 4.1,  0.70, 1.6e-2, 1.1, 4.2e-5, 2.2], \
            [180, 1.1e12, 10,  0.5, 4.0,  0.71, 8.7e-3, 1.2, 1.6e-5, 2.2], \
            [180, 1.1e12, 10,  0.3, 3.4,  0.71, 1.1e-2, 1.2, 2.1e-5, 2.2], \
            [150, 6.3e11, 10,  0.5, 6.1,  0.68, 2.8e-3, 1.2, 7.4e-6, 2.3], \
            [150, 6.3e11, 10, 0.3, 4.9,  0.68, 3.4e-3, 1.1, 9.8e-6, 2.2], \
            [150, 6.3e11, 10, 0.1, 3.0, 0.69, 8.1e-3, 1.2, 1.9e-5, 2.2], \
            [120, 3.2e11, 10, 0.5, 6.1,  0.69, 1.4e-3, 1.2, 2.3e-6, 2.2], \
            [120, 3.2e11, 10, 0.3, 5.0,  0.70, 1.9e-3, 1.2, 2.9e-6, 2.2], \
            [120, 3.2e11, 10, 0.1, 3.1,  0.71, 3.7e-3, 1.2, 5.1e-6, 2.2], \
            [120, 3.2e11, 20, 0.5, 9.6,  0.69, 7.0e-4, 1.2, 1.2e-6, 2.2],\
            [120, 3.2e11, 20, 0.3, 7.9,  0.70, 9.5e-4, 1.2, 1.5e-6, 2.2],\
            [120, 3.2e11, 20, 0.1, 4.9,  0.71, 1.9e-3, 1.2, 2.7e-6, 2.2]
            ])

        # Parameterization of percipitation model used in CGMBrush.  
        # XRvir is how many virial radii to go out for extended profile (this makes it so integrates to total mass within this radius)
        self.XRvir = XRvir 

        # Zmetal is the metalicity
        # 0.1 is only for lower mass halos, to use it will need to modify the code.
        if (Z_METAL != 0.3 and Z_METAL != 0.5):
            raise Exception("Invalid Z_METAL value, choose 0.3 or 0.5.")
        self.Z_METAL = Z_METAL  

        # cooling time to dynamical time ratio that specifies model. This is only option for below table. See Voit et al 2018 for more options.
        self.TRATCRIT = 10  

        # chooses metalicity and cooling criterion parameters we will interpolate on
        self.reducedarr = self.fitarray[(self.fitarray[:, 3] == self.Z_METAL) & (self.fitarray[:, 2] ==  self.TRATCRIT)]
        self.reducedarr = self.reducedarr[::-1] #reverses array       


    #density profile in comoving units
    #Mvir is the halo mass
    def get_analytic_profile(self, Mvir: float, redshift: float):
        #rvals = np.logspace(-3, 2, num=200, base=10)

        comoving_rvir= halo.comoving_rvir(cosmo, Mvir, redshift)
        rvirkpc = KPCINMPC * comoving_rvir
        epsilon = 1  #kpc -- shouldn't trust model at smaller values so always soften with this scale 

        #parameters of our percipitation profile
        n1,n2,xi1,xi2,neconstant, XRvir, rmax = self.get_precipitation_params(np.log10(Mvir), rvirkpc, redshift)
        rcomoving_kpc = np.logspace(1, np.log10(XRvir*rvirkpc), num=500, base=10)
        #neconstant = 0
        #print("neconst = ", neconstant, n1,n2,xi1,xi2, XRvir, " rmax = ", rmax, rvirkpc)

        final_ar = self.precipitation_func_array(rcomoving_kpc/(1+redshift), n1,n2,xi1,xi2,neconstant, XRvir*rvirkpc/(1+redshift), rmax, epsilon)
        

        return rcomoving_kpc/KPCINMPC, final_ar/(1+redshift)**3  #care to go back to comoving quantities

    #profile in physical units at z=0 (from appendix of https://arxiv.org/pdf/1811.04976.pdf) in physical distance units (In contrast to all other parts of code)
    #the constant density out to XRvir times the virial radius is so that the total mass in baryons is included
    def precipitation_func_array(self, rphyskpc, n1, n2, xi1, xi2, neconstant, tophat_limit, rmax, epsilon):

        # Make sure the power law contributes only up to rmax, and tophat only up to XRvir * Rvir
        points_in_tophat_limit = (np.array(rphyskpc) <= tophat_limit).astype(int) # makes an array for each point in the rphyskpc array, 1 if inside limit, 0 for outside
        points_in_powerlaw_limit = (np.array(rphyskpc) <= rmax).astype(int) # similar to above but for the rmax cutoff for power law part of profile
        
        return points_in_powerlaw_limit/np.sqrt((n1*(rphyskpc+epsilon)**-xi1)**-2 + (n2*((rphyskpc+epsilon)/100)**-xi2)**-2) + points_in_tophat_limit*neconstant

    #profile in physical units at z=0 (from appendix of https://arxiv.org/pdf/1811.04976.pdf) in physical distance units (In contrast to all other parts of code)
    #the constant density out to XRvir times the virial radius is so that the total mass in baryons is included
    def precipitation_func(self, rphyskpc, n1, n2, xi1, xi2, neconstant, rmax, epsilon):

        result = neconstant
        if rphyskpc <= rmax:
            result += 1/np.sqrt((n1*(rphyskpc+epsilon)**-xi1)**-2 + (n2*((rphyskpc+epsilon)/100)**-xi2)**-2)
        return result

    # Outputs percipitation model parmameters plus the constant density for the tophat
    def get_precipitation_params(self, log10Mhalo: float, comoving_rvir_kpc: float, redshift: float, calc_neconst_flag = True):
        Mvir = 10**log10Mhalo
        M_bary = Mvir * cosmo.fb
        log10Mhalo_z0 = log10Mhalo + 3/2*np.log10(1+redshift)  #This is how the Voit profile maps in redshift (the gas profile is fixed at vcir)
        
        rvir_physkpc = comoving_rvir_kpc/(1+redshift)

        # Interpolation of the Table in Voit paper
        logn1 = interp1d(np.log10(self.reducedarr[:, 1]), np.log10(self.reducedarr[:, 6]), kind='linear', fill_value='extrapolate')(log10Mhalo_z0)   
        n1 = 10**logn1
        xi1 = interp1d(np.log10(self.reducedarr[:, 1]), self.reducedarr[:, 7], kind='linear', fill_value='extrapolate')(log10Mhalo_z0) 
        logn2 = interp1d(np.log10(self.reducedarr[:, 1]), np.log10(self.reducedarr[:, 8]), kind='linear', fill_value='extrapolate')(log10Mhalo_z0)   
        n2 = 10**logn2
        xi2 = interp1d(np.log10(self.reducedarr[:, 1]), self.reducedarr[:, 9], kind='linear', fill_value='extrapolate')(log10Mhalo_z0) 

        #Calculates the constant density need to conserve mass assuming this extends to XRvir times the virial radius
        neconst = 0
        if calc_neconst_flag == True:
            r_physkpc = np.logspace(0, np.log10(self.XRvir*rvir_physkpc), 500)#radial bin array 
            #print("max r considered = ",  r_physkpc[-1], rvir_physkpc)
            #Voit 2018 fitting formulae 
            rhoarr = np.array(1/np.sqrt((n1*(r_physkpc)**-xi1)**-2 + (n2*(r_physkpc/100)**-xi2)**-2))

            #Integrate to see how much mass is missed by this profile  (I've checked these seems reasonable)
            #     rhointerp = interp1d(np.log(rhoarr[0]), 4.*np.pi*rhoarr[0]**3*rhoarr[1], kind='cubic', fill_value='extrapolate')
            # Logarithmic integral: 4 pi r^3 dlog(r) = 4 pi r^2 dr
            rhointerp = interp1d(np.log(r_physkpc), 4*np.pi*r_physkpc**3*rhoarr, kind='cubic') 

            #print("n1 n2, xi1, xi2", n1, n2, xi1, xi2)
            #print("rhoarr at 10 kpc = ", rhointerp(np.log(10))/( 4.*np.pi*10**3))
            #print("rhoarr at 100 kpc = ", rhointerp(np.log(100))/( 4.*np.pi*10**6))
            
            conv = (msun/(mu*mprot))/kpc**3
            mtotal = 0
            rmax = self.XRvir*rvir_physkpc*1.1 #set to just a little larger
            
            for i in range(1, len(r_physkpc)):
                mtotal += integrate.quad(rhointerp, np.log(r_physkpc[i-1]), np.log(r_physkpc[i]))[0]/conv
                if mtotal > M_bary:
                    if self.debug:
                        print(" Total mass fulfilled early")
                    rmax = r_physkpc[i]
                    break
            #print("mtot percip in 1e12 = ", mtotal/cosmo.fb/1e12, 10**(log10Mhalo-12))
        
            # Add in rest of mass (this should compensate for actual halo mass and not rescaled)
            if rmax > self.XRvir*rvir_physkpc:
                neconst = (M_bary-mtotal) / (4.*np.pi/3.*(self.XRvir*rvir_physkpc)**3) * conv    
        
            if self.debug:
                print(" Mass fraction in non-tophat: {:.3f}".format(mtotal/M_bary))

        return n1, n2, xi1, xi2, neconst, self.XRvir, rmax
        #return rkpc, rhoarr

    #outputs mask; All length scales have to be converted into units of cellsize
    def get_mask(self, mass: float, comoving_rvir_Mpc: float, redshift: float, resolution: int, cellsize: float, fine_mask_len: int):

        y,x = np.ogrid[-1*fine_mask_len: fine_mask_len, -1*fine_mask_len: fine_mask_len] # shape is (1,40*res) and (40*res,1)

        comoving_rvir_kpc = KPCINMPC * comoving_rvir_Mpc
        cellsize_kpc = KPCINMPC *cellsize  # kpc

        #parameters of our percipitation profile
        n1,n2,xi1,xi2,neconstant, XRvir, rmax = self.get_precipitation_params(np.log10(mass), comoving_rvir_kpc, redshift, True)
        
        #     f1= lambda x, y, z: my_func(((x**2+y**2+z**2)**.5), n1,n2,xi1,xi2,neconstant,cellsize_kpc)

        #integrate to project to 2D

        epsilon = self.epsilon*cellsize_kpc
        virial_radius = comoving_rvir_kpc/(1+redshift)
        tophat_limit = XRvir*virial_radius
        func = lambda x, y, z: self.precipitation_func(np.sqrt(x**2+y**2+z**2)*cellsize_kpc/(1+redshift), n1, n2, xi1, xi2, neconstant, rmax, epsilon)
                            
        interation_max = tophat_limit / cellsize_kpc
        mask = project_spherical_3Dto2D_optimized(func,x,y,interation_max)
        return mask.astype(float)



# This function subtracts the halos from the density field

# arguments: 
# haloArray: dataframe of halos, sorted by mass.
# bin_markers: array giving the indexes of halos at the edges of mass bins
# profile: tophat, NFW etc
def subtract_halos(provider: SimulationProvider, haloArray, bin_markers, profile: CGMProfile, redshift: float, halo):
    
    # TODO I think this is effectively hardcoded to the 256 Bolshoi grid size.
    df = haloArray
    no_cells = provider.halofieldresolution # the resolution multiplier isn't used for subraction, only addition
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

    coarse_cellsize = provider.Lbox/provider.halofieldresolution 
    fine_cellsize = coarse_cellsize / scale_down 
    
    # loops through the list of dataframes each ordered by ascending mass
    for j in range(0,chunks):

        if bin_markers[j] == bin_markers[j+1]:
            # Special case - mass bin is empty. Just set results to 0 for this mass bin.
            Mvir_avg[j] = 0
            conv_rad[j] = 0
            convolution[j,:,:] = 0
            continue

        Mvir_avg[j] = np.mean((df['Mvir'][bin_markers[j]:bin_markers[j+1]]))/cosmo.h
        conv_rad[j] = halo.comoving_rvir(cosmo, Mvir_avg[j], redshift) # comoving radius

        fine_mask = profile.get_mask(Mvir_avg[j], conv_rad[j], redshift, 1, fine_cellsize, fine_mask_len)

        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
        
        # Area of cells needed for normalization
        totalcellArea4 = sum(sum(coarse_mask))* ((coarse_cellsize)**2)

        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])   
        
        # Translate coordinates for comomving Mpc to halo grid coordinates.
        ix = ((np.rint(no_cells*((df[bin_markers[j]:bin_markers[j+1]]['x'].values)/(provider.Lbox*cosmo.h))))%no_cells).astype(int)
        iy = ((np.rint(no_cells*((df[bin_markers[j]:bin_markers[j+1]]['y'].values)/(provider.Lbox*cosmo.h))))%no_cells).astype(int)
        xy=(ix,iy)

        # BUG issue: the method does not add repeated coordinates
        halo_cell_pos[xy] += 1
        
        # convolve the mask and the halo positions
        convolution[j,:,:] = (Mvir_avg[j]/(totalcellArea4)) * my_convolve(halo_cell_pos,coarse_mask)    
        
    return (convolution.sum(0))*(Mpc**-3 *10**6)*nPS*cosmo.fb

def T_vir(M, z):    
    """
    Calculates the virial temperature of a halo.

    Arguments:
    M: halo mass (solar masses)
    z: redshift
    """
    # TODO look into this contribution to virial temperature
    #k = 0.01720209895 #rad
    #H0 = 100 #h km s^−1 Mpc^−1
    #Omegak = -k/(H0**2)
    #OmegaMz = (OmegaM*(1+z)**3)/((OmegaM*(1+z)**3)+OmegaL+(Omegak*(1+z)**2))
    mu = 0.6
    
    temperature = 1.98e4 * (mu / 0.6) * pow(M / 1e8, 2/3) * ((1+z)/10) # * (OmegaM/OmegaMz)**(1/3)
    return temperature

def T_anisotropy(DM, T_vir, z):
    """
    Calculates the temperature anisotropy of a halo.
    
    Arguments:
    DM: dispersion measure (parsec cm^-3)
    T_vir: virial temperature (K)
    z: redshift
    """
    # CMB temperature accounting for redshift
    T_cmb = 2.725 * (1+z) #K
    e_charge = 4.8032e-10 #cm^3/2 g^1/2 s^-1
    m_e = 9.10938356e-28 #g
    c = 2.99792458e10 #cm s^-1
    k_b = 1.3807e-16 #cm^2 g s^-2 K^-1
    
    #sigma_T = 8*np.pi/3 * pow(e_charge**2/(m_e*c**2), 2)
    sigma_T = 6.6524587158e-25 #cm^2
    
    dT = DM * kpc/10**3 * T_cmb * T_vir * sigma_T * 2 * k_b / (m_e * c**2)
    
    return dT


def normalize_and_convert_field(field, mass, mass_area):
    """Normalizes a field so the hydrogen fraction adds up to the mass provided. Converts to [pc cm^-3] units."""
    return (mass/mass_area) * (Mpc**-3 * 10**6) * nPS * cosmo.fb * field


def normalize_and_convert_mask(mask, mass, cellsize):
    """Normalizes a mask so the hydrogen fraction adds up to the mass provided. Converts to [pc cm^-3] units.
    Returns the mask and the total mass of the summed, unnormalized mask."""
    mass_area = sum(sum(mask)) * ((cellsize)**2)
    if mass_area != 0:
        return normalize_and_convert_field(mask, mass, mass_area), mass_area
    return mask, mass_area


def convolve_DM_for_bin(halo_cell_pos, mask, cellsize, Mvir_avg, redshift):
    """
    Convolves a 2D field of halo positions and a 2D mask of weights for a given mass bin.
    """
    
    # store addition masks
    final_mask, mass_area = normalize_and_convert_mask(mask, Mvir_avg, cellsize)
        
    if mass_area != 0:
        # Normalize and weight the results appropriately. 
        convolved_field = normalize_and_convert_field(my_convolve(halo_cell_pos, mask), Mvir_avg, mass_area)

        return convolved_field, final_mask
    else:
        # Empty mask means mass bin is empty. Skip convolution and return 0's
        return np.zeros(halo_cell_pos.shape), mask 

def convolve_dT_for_bin(halo_cell_pos, mask, cellsize, Mvir_avg, redshift):

    # TODO this hasn't been checked in a while
    Tvir_avg = T_vir(Mvir_avg, redshift)
    totalcellArea4 = sum(sum(mask)) * ((cellsize)**2)
    
    if totalcellArea4 != 0:
        # convolve the mask and the halo positions
        convolution = T_anisotropy((Mvir_avg/totalcellArea4) * (Mpc**-3 *10**6) * nPS * cosmo.fb, Tvir_avg, redshift) * my_convolve(halo_cell_pos, mask)
                    
        # store addition masks
        final_mask = T_anisotropy((Mvir_avg/totalcellArea4) * (Mpc**-3 *10**6)*nPS*cosmo.fb, Tvir_avg, redshift) * mask
        
        return convolution, final_mask
    else:
        # If the mass bin is empty, then skip convolution and return 0's
        return np.zeros(halo_cell_pos.shape), mask


def make_halo_DM_map(provider: SimulationProvider, haloArray, resolution: int, bin_markers, profile: CGMProfile, redshift: float):
    """
    Creates a map of dispersion measure.

    Arguments:
    haloArray: DataFrame of halos, sorted by mass.
    resolution: resolution of the halo field
    bin_markers: array giving the indexes of halos at the edges of mass bins
    profile: tophat, NFW etc
    redshift: redshift of halos
    """

    convolved_field, conv_rad, addition_masks, Mvir_avg = add_halos(provider, haloArray, resolution, bin_markers, profile, redshift, convolve_DM_for_bin)
            
    return convolved_field, conv_rad, addition_masks, Mvir_avg

def make_halo_dT_map(provider: SimulationProvider, haloArray, resolution: int, bin_markers, profile: CGMProfile, redshift: float):
    """
    Creates a map of CBM temperature anisotropy from the provided halos.
    
    Arguments:
    haloArray: dataframe of halos, sorted by mass.
    resolution: resolution of the halo field
    bin_markers: array giving the indexes of halos at the edges of mass bins
    profile: tophat, NFW etc
    redshift: redshift of halos
    """

    convolution, conv_rad, addition_masks, Mvir_avg = add_halos(provider, haloArray, resolution, bin_markers, profile, redshift, convolve_dT_for_bin)
    Tvir_avg = T_vir(Mvir_avg, redshift) # TODO duplicate calculation 

    return convolution, conv_rad, addition_masks, Tvir_avg


def add_halos(provider: SimulationProvider, haloArray, resolution: int, bin_markers, profile: CGMProfile, redshift: float, per_bin_func):
    """
    Performs a convolution between halo positions (via the haloArray parameter) and profile (via the profile parameter).

    Arguments:
    haloArray: DataFrame of halos, sorted by mass.
    resolution: int of resolution, will be multiplied by 1024
    bin_markers: array giving the indexes of halos at the edges of mass bins
    profile: tophat, NFW etc
    """
    
    no_cells = provider.halofieldresolution * resolution
    bins = len(bin_markers) - 1

    # array of halo masses and radii
    Mvir_avg = np.zeros(bins)
    conv_rad = np.zeros(bins)

    # convolution results
    convolution_summed = np.zeros([no_cells,no_cells])
    
    # fine and coarse map settings
    # TODO cleanup this
    fine_mask_len = 20*resolution # TODO CGMBrush is forcing this mask size scheme. Should implementers get to choose?
    scale_down = 2  # making the grid coarser    
    nbig = fine_mask_len*2 # fine masks are 40*resolution by 40*resolution
    nsmall = nbig//scale_down # coarse masks are half wide and long

    coarse_cellsize = provider.Lbox / no_cells
    fine_cellsize = coarse_cellsize / scale_down

    # store all (coarse) profile masks
    addition_masks =[]

    # loops through the list of dataframes each ordered by ascending mass
    for j in range(0,bins):

        if bin_markers[j] == bin_markers[j+1]:
            # Special case - mass bin is empty. Just set results to 0 for this mass bin.
            Mvir_avg[j] = 0
            conv_rad[j] = 0
            addition_masks.append(np.zeros((nsmall,nsmall))) # TODO make smaller, just doing this for safety for now
            continue

        Mvir_avg[j] = np.mean((haloArray['Mvir'][bin_markers[j]:bin_markers[j+1]])) / cosmo.h  #Matt: Would be much better to put in h at time we read in file
        conv_rad[j] = halo.comoving_rvir(cosmo, Mvir_avg[j], redshift) # comoving radius

        fine_mask = profile.get_mask(Mvir_avg[j], conv_rad[j], redshift, resolution, fine_cellsize, fine_mask_len)

        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)

        # trim away 0s around mask to make it as small as possible (speeds up convolution)
        with np.printoptions(precision=3, linewidth=1000, threshold=sys.maxsize):
            #print(coarse_mask.shape)
            i = 0
            length = coarse_mask.shape[0]
            assert (coarse_mask.shape[0] == coarse_mask.shape[1])

            # test if edges of mask are all zeros
            while not (np.any(coarse_mask[i]) or np.any(coarse_mask[length-1-i]) or np.any(coarse_mask[:,0]) or np.any(coarse_mask[:,i])):
                i+=1
                
            if i > 0:
                to_del = np.concatenate((np.arange(0,i), np.arange(length-i,length)))
                #print(to_del)
                coarse_mask = np.delete(np.delete(coarse_mask, to_del, axis=1), to_del, axis=0) 

            #print(coarse_mask)
                    

        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])    
        
        # Translate the x,y values from comoving coordinates to fine grid coordinates
        # The modulus just ensures that halos on the edge that are rounded up wrap around to the 0th index location
        ix = ((np.rint(no_cells*((haloArray[bin_markers[j]:bin_markers[j+1]]['x'].values)/(provider.Lbox*cosmo.h))))%no_cells).astype(int)
        iy = ((np.rint(no_cells*((haloArray[bin_markers[j]:bin_markers[j+1]]['y'].values)/(provider.Lbox*cosmo.h))))%no_cells).astype(int)
        xy=(ix,iy)

        # Adnan: BUG issue: the method does not add repeated coordinates. Ian: Is that right?
        halo_cell_pos[xy] += 1

        # convolve the mask and the halo positions
        convolution, final_mask = per_bin_func(halo_cell_pos, coarse_mask, coarse_cellsize, Mvir_avg[j], redshift) 
        
        convolution_summed += convolution

        # store addition masks
        addition_masks.append(final_mask)
        
    return convolution_summed, conv_rad, addition_masks, Mvir_avg

# Halos removed field

#This function combines many steps of subtraction and smoothing to yield a density field from which halos have been removed
def halos_removed_field(provider: SimulationProvider, current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,subtraction_halo_profile,resolution,sigma_gauss,width_sinc, halo):
    
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    bin_markers= halo_array_for_convolution[1]
    
    # convolve halos
    subtraction_profile = subtract_halos(provider, df,bin_markers,subtraction_halo_profile,redshift, halo)
    assert not np.any(np.isnan(subtraction_profile))
    subtraction_profile_smooth = gauss_sinc_smoothing(subtraction_profile,sigma_gauss,width_sinc,1, provider.halofieldresolution)
    assert not np.any(np.isnan(subtraction_profile_smooth))
    
    # create coarse grid
    subtraction_coarse= smoothfield(subtraction_profile_smooth, provider.halofieldresolution, den_grid_size)
    assert not np.any(np.isnan(subtraction_coarse))
    
    # remove halos from the density field
    halos_removed_coarse=removeConvolvedHalos(density_field,subtraction_coarse)
    assert not np.any(np.isnan(halos_removed_coarse))
    
    return halos_removed_coarse,subtraction_coarse,subtraction_profile,subtraction_profile_smooth


# Function subtracts and adds halos
def convolution_all_steps_final(provider, current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,halos_removed_coarse,
                       addition_halo_profile: CGMProfile,resolution,sigma_gauss,width_sinc):
    
    # setup inputs for convolution
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    bin_markers= halo_array_for_convolution[1]
    
    # convolve halos for adding back
    addition_profile_initial=make_halo_DM_map(provider, df,resolution,bin_markers,addition_halo_profile,redshift)
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
    

def halo_subtraction_addition(sim_provider : SimulationProvider,den_grid_size,RS_array,min_mass,max_mass,log_bins,subtraction_halo_profile: CGMProfile,
                             addition_profile: CGMProfile,resolution, halo):
    """
    This function runs the convolution code (subtraction and addition) and returns all the results. It will run for every redshift given serially.

    Outputs: a tuple with halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, halos removed field, virial radii, halo masses, subtraction field, subtraction field smoothed
    """

    # Details of halo profiles

    # sigma is the variance of the gaussian. sigma**2 = 16 or 3*sigma**2 = 16
    sigma_gauss = 4

    # nsinc is the full width of the box function. So 4 means two boxes on either side of the point.
    width_sinc = 4

    grid_length = sim_provider.halofieldresolution * resolution

    halos_reAdded = np.zeros([len(RS_array),grid_length,grid_length])
    
    halos_subtraction_coarse = np.zeros([len(RS_array),den_grid_size,den_grid_size])
    halo_field = np.zeros([len(RS_array),grid_length,grid_length])
    halo_masks = [] # list (for each z) of lists of 2D arrays of masks (ragged sizes)
    halos_removed_fields = np.zeros([len(RS_array),den_grid_size,den_grid_size]) # This should be the same as fine_mask_len in add_halos function
    
    for i in range(0, len(RS_array)):
        
        redshift = RS_array[i]
        density_field = sim_provider.get_density_field(redshift, den_grid_size)
        halos = sim_provider.get_halos(redshift)

        halos_removed = halos_removed_field(sim_provider, halos,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,subtraction_halo_profile,resolution,sigma_gauss,width_sinc, halo)
              
        conv_all_steps = convolution_all_steps_final(sim_provider, halos,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,halos_removed[0],
                           addition_profile,resolution,sigma_gauss*resolution,width_sinc)
        
        
        halos_reAdded[i,:,:] = conv_all_steps[0]
        
        halos_subtraction_coarse[i,:,:] = halos_removed[1]
        halo_field[i,:,:] = conv_all_steps[2]
        halo_masks.append(conv_all_steps[3])
        halos_removed_fields[i,:,:] =  halos_removed[0]
        
        
    # returns halo array and masks used to add halos back
    return halos_reAdded,halo_masks,halos_subtraction_coarse,halo_field,halos_removed_fields, conv_all_steps[5],conv_all_steps[6],halos_removed[2],halos_removed[3]



####################################
# DM profile for a mass bin
####################################

# Create radial profile of a field
def radial_profile(data, center_x,center_y):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    #print(r)
    r = r.astype(np.int)
    #print(r)
    #exit()
    
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 


def make_halo_square(DM_field, ix, iy, crop_grid):
    """Creates a square cutout of the DM field from around the ix, iy point 
    with dimensions crop_grid x crop_grid."""

    res=DM_field.shape[0]    
    trimmed = -1 # sentinal value to caller to drop this halo square
     
    # TODO This ignores halos near the edge; we can wrap around and not drop these halos.
    if ix > int(crop_grid) and ix < res - int(crop_grid) and iy > int(crop_grid) and iy < res- int(crop_grid):
        trimmed = DM_field[ix-int(crop_grid/2):ix+int(crop_grid/2),iy-int(crop_grid/2):iy+int(crop_grid/2)]

    return trimmed


# Single function for radial profile of DM for a given DM array and grid size
def DM_vs_radius(DM_field, halo_data_frame, crop_grid_dim, bin_markers, provider: SimulationProvider):    
    """Creates radial profile of the DM for halos in the DM_field provided. Uses the halo_data_frame, which is a DataFrame and must be sorted in ascending mass order,
     alongside the bin_markers which provides the indexes where mass bins change in the halo_data_frame. crop_grid_dim is used to choose how large a square around the halo
     centers to cutout. """

    # TODO: Matt says we can calculate using every 1/4 the pixels or something.
    
    assert crop_grid_dim % 2 == 0, "The present implementation requires crop_grid_dim to be even." # BUG We should probably make it odd actually.

    num_bins = len(bin_markers) -1 
    center = int(crop_grid_dim/2)
    
    DM_mass_bin= np.zeros([num_bins, crop_grid_dim, crop_grid_dim])
    rad_prof_mass_bin = []

    no_cells=DM_field.shape[0]
    # Translate coordinates from comoving Mpc to fine grid coordinates
    ix = ((np.rint(no_cells*((halo_data_frame['x'].values)/(provider.Lbox*cosmo.h))))%no_cells).astype(int)
    iy = ((np.rint(no_cells*((halo_data_frame['y'].values)/(provider.Lbox*cosmo.h))))%no_cells).astype(int)

    #print("ix: {}".format(ix))
    #print("iy: {}".format(iy))

    # Go through each mass bin
    for i in range(0, num_bins):
        # get a square cutout from the DM field around the center of each halo within this mass bin
        num_halos = bin_markers[i+1] - bin_markers[i]
        #print("Mass bin {}: {} halos".format(i,num_halos))

        for halo_index in range(bin_markers[i], bin_markers[i+1]):
            halo_square = make_halo_square(DM_field, ix[halo_index], iy[halo_index], crop_grid_dim)
            #print(halo_square)
            if isinstance(halo_square, (list, tuple, np.ndarray)):
                DM_mass_bin[i,:,:] = DM_mass_bin[i,:,:] + halo_square
            else:
                # currently we drop halos near the edge that we don't have room to cutout the halo square
                # TODO if easy, use periodic boudnaries to calculate it all correctly.
                num_halos -= 1 
           
        #print("Bin {}: {} usable halos.".format(i, num_halos))
             
        if (num_halos <= 0):
            num_halos = 1

        DM_mass_bin[i,:,:] = DM_mass_bin[i,:,:] / num_halos # finishing computing average
        rad_prof_mass_bin.append(radial_profile(DM_mass_bin[i,:,:],center,center))
    
    arr = np.array(rad_prof_mass_bin)
    return arr, DM_mass_bin
    
def profile_of_masks(mask_array):
    """This calculate a radial profile of all the masks provided in mask_array, seperately processing each mass bin into an average values by radius."""

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

    print("saving output to ", file_path)
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
    #print("Trying to load ", file_path)
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
    def __init__(self, addition_profile: CGMProfile, provider: SimulationProvider = None, folder=VAR_DIR, resolution=1, den_grid_size=256, RS_array=[0], datestamp=str(datetime.date.today())):
        
        # Profile to use for adding in CGM
        self.addition_profile = addition_profile
        self.provider = provider
        self.resolution = resolution # actual fine grids are generally this x provider.halofieldresolution
        self.folder = folder

        # Resolution: choose between 256 and 512 grid TODO this is Bolshoi specific
        if den_grid_size != 256 and den_grid_size != 512:
            raise ValueError("Only resolutions 256 and 512 are allowed")
        self.den_grid_size = den_grid_size 

        # User provides a redshift array
        self.RS_array = RS_array # For a single box, we only use the redshift 0 box

        # Profile used for subtracting halos from the density field
        self.subtraction_halo_profile = NFWProfile()

        # TODO these values change the results but are not included in the filename for saved data,
        # so two configuration that have different values here but are otherwise the same will collide.
        # TODO: best fix is probably to make saved files have a hashed value for all parameters of the
        # config that affect them. Fold in all the previo stuff into that. CGMProfiles must be hashable.
        self.min_mass = DEFAULT_MIN_MASS # halos smaller than this shouldn't have much of a CGM
        self.max_mass = DEFAULT_MAX_MASS # this is a little bigger than the biggest halo in Bolshoi
        self.log_bins = DEFAULT_MASS_BIN_COUNT + 1 # this means 60 bins; TODO make this more intuitive... 
        self.datestamp = datestamp
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
        """
        Constructs a filename for this configuration. It is important that this precription stays stable. 
        
        Unfortunately the present implementation does not account for every parameter and so this is imperfect, but works for most purposes.
        """
        scaling = ''
        if self.addition_profile.parameter_name_str: 
            scaling = '_' + str(self.addition_profile.parameter_name_str)
        z_str = ''
        if self.RS_array != [0]:
            z_str = '_z'
            for z in self.RS_array: # z_0.1_0.2_0.3 for example
                z_name = self.provider.get_z_name(z)
                z_str += '_' + z_name
        return self.addition_profile.name + str(self.resolution) + scaling + '_' + str(self.den_grid_size) + z_str + "_" + self.datestamp

    def convert_and_save(self):
        """Makes results available as a dictionary, which is how the .npz archive is formatted. 

        halo_subtraction_addition returns a tuple. Reading .npy files gets a tuple or array. We want everything converted to dictionary and npz.
        
        Note that .npz archives are loaded lazilly, which is important for managing memory usage."""
        # halo_subtraction_addition returns a tuple. Reading .npy files gets a tuple or array. 
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

    def get_add_mask_for_mass(self, mass, z_index):
        """Selects the addition mask closest to the specified halo mass (in solar masses) for the data from the given redshift's index and returns it."""
        masses = self.get_halo_masses()
        index = np.argmin(np.abs(masses - mass))
        masks = self.get_addition_masks()
        return masks[z_index][index]

    def run(self, trace=False, results_in_memory=True, load_from_files=False):
        """Run convolutions for this configuration."""

        filename = self.get_filename()
        file_path = os.path.join(self.folder, filename + ".npz")

        #print("load_from_files = ", load_from_files)
        if load_from_files:
            try:
                print("Trying to load data from ", file_path, end=" ")
                self.npz = np.load(file_path, allow_pickle=True)
                print("done.")
                            
            except IOError:
                print("Previous results '" + filename + "' not found. Calculations will be made.")
                #pass # previous results file for this config doesn't exist, swallow and compute it instead

        if self.npz is None:
            print("Performing calculations for {}... ".format(filename), end="")
                                   
            if trace:
                pr = cProfile.Profile()
                pr.enable()

            self.results = halo_subtraction_addition(self.provider, self.den_grid_size, self.RS_array, self.min_mass, 
                                                self.max_mass, self.log_bins, self.subtraction_halo_profile, 
                                                self.addition_profile, self.resolution, halo)
            
            self.convert_and_save()
            self.npz = np.load(file_path, allow_pickle=True)
            print("done.")

            if trace:
                pr.disable()
                s = io.StringIO()
                ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
                ps.print_stats()

                version = 5
                # TODO auto version incrementing
                perf_file_path = os.path.join(VAR_DIR, 'perf_convo_' + filename + '_v%s.txt' % version)
                with open(perf_file_path, 'w') as f:
                    f.write(s.getvalue())
        
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

        #translated_file = self.get_filename() + "_translated"
        stacked_file = self.get_filename() + "_stacked"
        file_path = os.path.join(self.folder, stacked_file + ".npz")
    
        if load_from_files:
            try:
                print("Loading stacked fields... ", end="")
                self.stacked_npz = np.load(file_path, allow_pickle=True)
                print("done")        
            except IOError:
                print("Previous results '" + stacked_file + "' not found. Calculations will be made.")

        if self.stacked_npz is None:

            if len(self.RS_array) > 1:
                print("Creating Stacked Fields... ", end="")

                all_orig_fields = np.zeros((len(self.RS_array), self.den_grid_size, self.den_grid_size))
                for i in range(0, len(self.RS_array)):
                    all_orig_fields[i] = self.provider.get_density_field(self.RS_array[i], self.den_grid_size)
                translated_field = translate_field_stack(all_orig_fields, self.RS_array, self.seed)
                self.stacked_orig_field = np.zeros(translated_field.shape)
                for i in range(0, len(self.RS_array)):
                    self.stacked_orig_field[i,:,:] = redshifted_DM(translated_field[i,:,:], self.RS_array[i])

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
                self.stacked_npz = np.load(file_path, allow_pickle=True)

                print("done")
            else:
                raise ValueError('Generating a stacked field is only applicable with data from multiple redshifts.')

        if not results_in_memory:
            self.stacked_orig_field = None
            self.stacked_removed_field = None
            self.stacked_addition_field = None
            self.stacked_final_field = None
        

    def generate_DM_vs_radius_profile(self, load_from_files=False, index=0):
        """Generates DM vs Radius profiles. By default this uses the 1st box available, other boxes may be chosen by specifying the index.

        Note that the DM_vs_R1 variable only holds the most recently generated (or loaded) DM vs Radius info. So if two calls are made to this method
        for different boxes, the DM_vs_R1 field will contain the results of the second call only. It is better to use the data returned from this method.
        """
        self.DM_vs_R1 = None

        profile_file = '%s_DMvsR_prof' % self.get_filename()
        if index > 0:
            profile_file = profile_file + '_for_box_{}'.format(index)

        if load_from_files:
            try:
                self.DM_vs_R1 = loadArray(profile_file, folder=self.folder)            
            except IOError:
                print("Previous results '" + profile_file + "' not found. Calculations will be made.")

        if self.DM_vs_R1 is None:       
            print("Generating DM vs R profile for box {}".format(index))
            df = create_halo_array_for_convolution(self.provider.get_halos(self.RS_array[index]), self.min_mass, self.max_mass, self.log_bins)

            trim_dim = int(10*self.resolution)

            self.DM_vs_R1 = DM_vs_radius(self.get_final_field()[index], df[0], trim_dim, df[1], self.provider) [0]

            saveArray(profile_file, self.DM_vs_R1, folder=self.folder)
  
        return self.DM_vs_R1
    
    def generate_profile_of_masks(self, load_from_files=False, index=0):
        """Generates mask profiles. By default this uses the 1st box available, other boxes may be chosen by specifying the index.

        Note that the mask_profiles variable only holds the most recently generated (or loaded) DM vs Radius info. So if two calls are made to this method
        for different boxes, the mask_profiles field will contain the results of the second call only. It is better to use the data returned from this method.
        """
        self.mask_profiles = None
        
        mask_file = '%s_masks' % self.get_filename()
        if index > 0:
            mask_file = mask_file + '_for_box_{}'.format(index)

        if load_from_files:
            try:
                self.mask_profiles = loadArray(mask_file, folder=self.folder)            
            except IOError:
                print("Previous results '" + mask_file + "' not found. Calculations will be made.")


        if self.mask_profiles is None:
            print("Generating Mask Profiles for box {}".format(index))

            # zoom in on middle of masks
            full_mask_len = self.resolution * 20
            zoom_start = full_mask_len // 4
            zoom_end = 3 * zoom_start
            self.mask_profiles = profile_of_masks(self.get_addition_masks()[index, :, zoom_start:zoom_end, zoom_start:zoom_end])
            saveArray('%s_masks' % self.get_filename(), self.mask_profiles, folder=self.folder)
        
        return self.mask_profiles

    def clear_results(self):
        """Clears memory-intense results from memory. Saved files are preserved. Results can be recovered quickly by running with load_from_files=True."""
        self.results = None
        self.final_field = None
        self.addition_field = None
        self.stacked_addition_field = None
        self.stacked_final_field = None
        self.stacked_removed_field = None
        self.stacked_orig_field = None
        # TODO is calling del(...) better? 
        #gc.collect()
        # should allow garbage collection to happen
        
