from __future__ import print_function 
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from numpy import core
import scipy.integrate as integrate
from numpy import genfromtxt
from random import randrange
import random
import time
from numpy.polynomial import Polynomial as P
import pandas as pd
import os
from scipy.ndimage.filters import convolve
import abc
import cProfile
import io
import pstats
import datetime
import math
from scipy.interpolate import interp1d

########################################
# Paramters, Constants, and Functions
########################################

# Constants
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

# TODO Bolshoi specific?
sims_folder= "../sims"
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

# BUG isn't it + 82, not - 82?
def rho_vir(z):
    return (18*pi**2 - 82*q(z) - 39*q(z)**2)*(rho_c*(OmegaL + OmegaM *(1+z)**3))


def Rvir_den(z):
    return (4/3 * np.pi * rho_vir(z))**(1/3) # physical, units in 1/r**3

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


class BolshoiProvider(SimulationProvider):
    """SimulationProvider implementation for Bolshoi simulations."""

    # TODO pass in filename in constructor?
    def __init__(self):
        # Associative array of (redshift, resolution) => 3D numpy grid 
        # This works very well with np.savez, which is much faster to read
        # than the Bolshoi density files
        self.density_fields = {}  
        self.halos = {}

    def get_density_field(self, redshift: float, resolution: int):
        """Gets the density field for the given redshift and grid resolution."""

        # If in memory, use it
        if (redshift, resolution) in self.density_fields:
            return self.density_fields[(redshift, resolution)]

        # If not, try fast reading from saved off numpy array            
        filename = "bol_den_field_" + str(redshift) + '_' + str(resolution) 
        result = []
        try:
            result = loadArray(filename)
        except IOError:
            pass 
        
        # Do the slow reading from the bolshoi file and save off fast copy
        if len(result) == 0:
            result = self.import_density_field(redshift, resolution)
            if len(result) == 0:
                raise IOError("There was a problem importing the Bolshoi density field for (" + str(redshift) + ", " + str(resolution) + ")")
            saveArray(filename, result)

        self.density_fields[(redshift, resolution)] = result
        return result

    def get_halos(self, redshift: float) -> pd.DataFrame:
        """Gets halo information for the given redshift."""

        # If in memory, use it
        if redshift in self.halos:
            return self.halos[redshift]

        # If not, try fast reading from saved off numpy array            
        filename = "bol_halos_" + str(redshift)
        result = []
        try:
            result = loadArray(filename)
            # The rest of the code expects this to be a DataFrame, so convert it back
            result = pd.DataFrame(result, columns = ['row_id','x','y','z','Mvir','Mtot','Rvir','ix','iy','iz'])
        except IOError:
            pass 
        
        # Do the slow reading from the bolshoi file and save off fast copy
        if len(result) == 0:
            result = self.extract_halos(redshift)
            if len(result) == 0:
                raise IOError("There was a problem importing the Bolshoi halos for redshift" + str(redshift))
            saveArray(filename, result)

        self.halos[redshift] = result
        return result

    # Function to import density and halo tables for a given redshift
    def import_density_field(self, redshift, resolution):
        
        if resolution != 256 and resolution != 512:
            raise ValueError("Only resolution 256 or 512 is supported.")

        resStr = str(resolution)

        # Have to hardcode z=0 table because of the unique column names
        if redshift == 0:
            # reading density field and halos data
            file_path= os.path.join(sims_folder, 'dens'+resStr+'-z-0.csv.gz')
            pdDens=pd.read_csv(file_path)

            # extracting columns
            pdDensN=pdDens[['Bolshoi__Dens'+resStr+'_z0__ix','Bolshoi__Dens'+resStr+'_z0__iy','Bolshoi__Dens'+resStr+'_z0__iz','Bolshoi__Dens'+resStr+'_z0__dens']]

            # 3D density array
            pdDensN=pdDensN.sort_values(['Bolshoi__Dens'+resStr+'_z0__ix','Bolshoi__Dens'+resStr+'_z0__iy','Bolshoi__Dens'+resStr+'_z0__iz'])
            tden = pdDensN['Bolshoi__Dens'+resStr+'_z0__dens'].values
            tden2 = np.reshape(tden,(resolution,resolution,resolution))
            return ((tden2+1).sum(2))*10**6* dx* elecD(z) /(1+z)**2

        else:
            name = 'dens'+resStr+'-z-0.'+str(redshift)+'.csv.gz'
            den = pd.read_csv(name)
            den2=den[['Bolshoi__Dens'+resStr+'__ix','Bolshoi__Dens'+resStr+'__iy','Bolshoi__Dens'+resStr+'__iz','Bolshoi__Dens'+resStr+'__dens']]

            # 3D density array
            den_sorted=den2.sort_values(['Bolshoi__Dens'+resStr+'__ix','Bolshoi__Dens'+resStr+'__iy','Bolshoi__Dens'+resStr+'__iz'])
            den_vals = den_sorted['Bolshoi__Dens'+resStr+'__dens'].values
            den = np.reshape(den_vals,(resolution,resolution,resolution))
            
            return normDM((den+1).sum(2),0)
            #return normDM((den+1).sum(random.randint(0,1)),0)

    # Create a single array of density fields for various redshifts
    # Density fields: Choose 256 or 512 field
    def extract_all_den_fields(self, RS_array, den_grid_size):

        all_den_fields = np.zeros([len(RS_array),den_grid_size,den_grid_size])

        for i in range(0, len(RS_array)):
            all_den_fields[i,:,:] = self.import_density_field(i, den_grid_size)
        
        return all_den_fields 
        
    # Extract halos for a given redshift
    def extract_halos(self, redshift):
        # TODO filenaming issues? Adnan had a mismatch and just called files 0.1, 0.2, etc
        name = 'halo-z-0.'+str(redshift)+'.csv.gz'    
        file_path= os.path.join(sims_folder, name)
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

    The length will be up to logchunks, but may be smaller due to removed bins (if a mass bin would be empty)."""
    halos = haloArray_minmax(pdHalos,M_min,M_max)
    df = halos.sort_values(by='Mvir',ascending=True)
    sorted_haloMasses = df['Mvir'].values
    histH, binsH = np.histogram(sorted_haloMasses,bins=np.logspace(np.log10(M_min),np.log10(M_max), logchunks))
    bins=np.append([0],histH.cumsum()).tolist()

    # Remove any empty bins
    i = 0
    while i < len(bins) - 1:
        if bins[i] == bins[i+1]:
            del bins[i+1]
        else:
            i += 1
    
    return df,bins


# The user can create their own function
# All length scales have to be converted into units of cellsize

def precipitation_func(r, n1,n2,xi1,xi2,neconstant,cellsize_kpc,Rvirkpc,XRvir,redshift):
        
    #final_ar = np.array([1/np.sqrt(1/(n1*((r+.5)*cellsize_kpc)**-xi1)**2 + 1/(n2*((r+.5)*cellsize_kpc/100)**-xi2)**2)])# + neconstant
    x = (np.array(r) <= XRvir*Rvirkpc)

    # TODO I got this line broken from Adnan and had to guess where the ) goes to make it valid
    final_ar =   np.array([1/np.sqrt(1/(n1*((r+.5)*cellsize_kpc/(1+redshift))**-xi1)**2 + 1/(n2*((r+.5)*cellsize_kpc/(100*(1+redshift))**-xi2)**2)) + x.astype(int)*neconstant])
       

    return final_ar



# nePercipitation(10)

# pass rkpc 



# Convert custom 3D function to 2D by integrating over z
def func3Dto2D(f,x,y,Rvir):
    return integrate.quad(f,-(Rvir**2-(x**2+y**2))**.5,(Rvir**2-(x**2+y**2))**.5,args=(x,y))[0]

# TODO dedupe
def fire3Dto2D(f,x,y,Rvir):
    return integrate.quad(f,-1*np.real((Rvir**2-(x**2+y**2))**.5),np.real((Rvir**2-(x**2+y**2))**.5),args=(x,y))[0] 
    

# Project NFW profile from 3D to 2D
def NFW2D(x,y,rho_nought,R_s,Rvir):
    offset=float(.1)
    return integrate.quad(lambda x, y, z: rho_nought/(((offset+(x**2+y**2+ z**2)**.5)/R_s)*(1+((x**2+y**2+ z**2)**.5)/R_s)**2),-1*np.real((Rvir**2-(x**2+y**2))**.5),np.real((Rvir**2-(x**2+y**2))**.5),args=(x,y))[0]

# Function creates a smaller grid from a larger grid by smoothing
def smoothfield(big, nbig,nsmall):
    
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


# This function subtracts the halos from the density field

# arguments: 
# haloArray: dataframe of halos, sorted by mass.
# resolution: 1024 or 2048
# bin_markers: array giving the indexes of halos at the edges of mass bins
# profile: tophat, NFW etc
# scaling_radius: scale radius for tophat halos
def subtract_halos(haloArray,resolution,bin_markers,profile,scaling_radius,redshift):
    
    df = haloArray
    no_cells = 1024
    cellsize = L/(1024) 
    chunks = len(bin_markers) - 1
    
    # array of halo masses and radii
    Mvir_avg = np.zeros(chunks)
    conv_rad = np.zeros(chunks)
    Rvir_avg = np.zeros(chunks)

    # convolution mask array
    convolution = np.zeros([chunks,no_cells,no_cells])    
    #new_conv = np.zeros([chunks,no_cells,no_cells])
    # creating a coarse map out of a fine mask
    
    # fine mask
    fine_mask_len = 10
    fine_lower= -1*fine_mask_len
    fine_upper= fine_mask_len
    y,x = np.ogrid[fine_lower: fine_upper, fine_lower:fine_upper]
    
    # coarse mask
    scale_down = 2  # making the grid coarser
    smooth_r=2
    
    coarse_mask_len = int(fine_mask_len/scale_down)
    fine_mask= np.zeros([2*coarse_mask_len,2*coarse_mask_len])
    
    # loops through the list of dataframes each ordered by ascending mass
    
    for j in range(0,chunks):
        Mvir_avg[j] = np.mean((df['Mvir'][bin_markers[j]:bin_markers[j+1]]))/h
        conv_rad[j] = (1+redshift)*((Mvir_avg[j])**(1/3) / Rvir_den(redshift)) # comoving radius
        
        # NFW 
        # add redshift to the function above
        R_s=0
        rho_nought=0
        
        R_s = conv_rad[j]/(halo_conc(redshift,Mvir_avg[j])*cellsize)  
        rho_nought = rho_0(redshift,Mvir_avg[j],R_s)
        
        
        # Why are we multiplying by scale_down
        if profile == 'tophat':
            r = (x**2+y**2)**.5
            fine_mask = r <= (scaling_radius*scale_down*conv_rad[j]/cellsize) 
            fine_mask=fine_mask.astype(float)

            
        # spherical tophat
        elif profile == 'tophat_spherical':
            r = (x**2+y**2)**.5
            fine_mask = r <= (scaling_radius*scale_down*conv_rad[j]/cellsize)
            
            fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
            Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
            fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
        
        elif profile == 'NFW':
            vec_integral = np.vectorize(NFW2D)
            fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j])

            r=(x**2+y**2)**.5 # * scale_down
            
            fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
            
            fine_mask=fine_mask.astype(float)
        
        
        elif profile == 'custom':
            
            # Functions for profile
            # Currently hardcoded but should allow the user to provide the function as inputs
            
            f1_1 = lambda x,y,z: 1/(x**2+y**2+z**2+.5)**.5
            f1_2 = lambda x,y,z: 1/(x**2+y**2+z**2+.5)
            
            # Radius of first profile:
            R= (scaling_radius*scale_down*conv_rad[j]/cellsize)/2 
            
            
            vec_integral = np.vectorize(func3Dto2D)
            
            mask1 = x**2+y**2 < (R/2)**2
            mask1 = mask1.astype(float)*vec_integral(f1_1,x,y,R/2)

            mask2_1 = x**2+y**2 >= (R/2)**2
            mask2_1= mask2_1.astype(float)
            mask2_2 = x**2+y**2 <= (R)**2
            mask2_2= mask2_2.astype(float)
            mask2=mask2_1*mask2_2
            mask2 = mask2*vec_integral(f1_2,x,y,R)


            fine_mask=mask1+mask2
            
        
        elif profile == 'custom_tophat':
            
            # Functions for profile
            # Currently hardcoded but should allow the user to provide the function as inputs
            
            # Radius of first profile:
            R= (scaling_radius*scale_down*conv_rad[j]/cellsize)
            
            
            f1_1 = lambda x,y,z: np.exp(-((x**2+y**2+z**2)/R**2)**30)
            
            vec_integral = np.vectorize(func3Dto2D)
            
            mask1 = vec_integral(f1_1,x,y,R)

            fine_mask=mask1
        
        
        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        nbig = fine_mask_len*2
        nsmall = int(nbig/scale_down)
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
        
        # Area of cells needed for normalization
        totalcellArea4=0
        totalcellArea4 = sum(sum(coarse_mask))* ((cellsize)**2)

        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])   
        
        # The coordinates are being multiplied by 4 to yield the halo coordinates on the 1024 grid
        ix = ((((np.around(4*((df[bin_markers[j]:bin_markers[j+1]]['x'].values)/(250/256))))))%(1024)).astype(int)
        iy = ((((np.around(4*((df[bin_markers[j]:bin_markers[j+1]]['y'].values)/(250/256))))))%(1024)).astype(int)
        
        xy=(ix,iy)

        # issue: the method does not add repeated coordinates
        halo_cell_pos[xy] += 1
        
        # convolve the mask and the halo positions
        c1 = convolve(halo_cell_pos,coarse_mask, mode='wrap')
        #c2 = my_convolve(halo_cell_pos,coarse_mask)
        convolution[j,:,:] = (Mvir_avg[j]/(totalcellArea4)) * c1

        # My fft convolution method gives results that are different from old method at a level > 1e-03 sometimes. TODO which is more accurate? Does it matter?
        #new_conv[j,:,:] = (Mvir_avg[j]/(totalcellArea4)) * c2
        #with np.printoptions(threshold=np.inf,linewidth=np.inf, precision=2):
        #    assert (np.allclose(c1, c2, rtol=1e-04, atol=1e-04)), "New convolution is not equivalent to old (in wrap mode) during subtract_halos for j=" + str(j) + ".\n" + str(coarse_mask)
        
        
    
    return (convolution.sum(0))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM), conv_rad, Rvir_avg, fine_mask, coarse_mask,halo_cell_pos


class CGMProfile(metaclass=abc.ABCMeta):
    """Interface used by cgmbrush that handles creating a CGM profile to convolve."""

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_mask') and
        callable(subclass.get_mask) or 
        NotImplemented)

    def get_fine_mask_len(self, resolution: int):
        """Gets the length of the fine mask this profile will use. TODO But actually it is double this..."""
        return 20*resolution

    @abc.abstractmethod
    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, *args):
        """Constructs and returns a 2D mask to convolve with the halo locations
           for the specified halo parameters and redshift."""
        raise NotImplementedError

class TophatProfile(CGMProfile):

    def get_mask(self, mass: float, comoving_radius: float, redshift: float, resolution: int, scaling_radius: int, cellsize: float, *args):

        # fine mask
        # fine mask size has to correspond to the size of the mask that I eventually trim
        fine_mask_len = self.get_fine_mask_len(resolution)
        fine_lower= -1*fine_mask_len
        fine_upper= fine_mask_len
        y,x = np.ogrid[fine_lower: fine_upper, fine_lower:fine_upper] # shape is (1,40*res) and (40*res,1)

        scale_down = 2  # making the grid coarser    
        
        # TODO Why are we multiplying by scale_down
        r = (x**2+y**2)**.5
        fine_mask = r <= (scaling_radius * scale_down * comoving_radius / cellsize) 
        return fine_mask.astype(float)

    
# From rhogas Fire, TODO reorganize
# The user can create their own function
# All length scales have to be converted into units of cellsize
def fire_func(r, rmax, Rinterp, rho0, cellsize):
    R1 = rmax/cellsize  # Convert length scale into units of cellsize
    R2 = Rinterp/cellsize
    return rho0*np.exp(-r/R1)*  ((r+.5)/R2)**-2

# Status: This is the latest convolution function

# The function add halos

# arguments: 
# haloArray: dataframe of halos, sorted by mass.
# resolution: int of resolution, will be multiplied by 1024
# bin_markers: array giving the indexes of halos at the edges of mass bins
# profile: tophat, NFW etc
# scaling_radius: scale radius for tophat halos
def add_halos(haloArray, resolution: int, bin_markers, profile: CGMProfile, scaling_radius: int, redshift: float):
    chunks = len(bin_markers) - 1
    df = haloArray
    no_cells = 1024 * resolution
    cellsize = L / no_cells
    
    # array of halo masses and radii
    Mvir_avg = np.zeros(chunks)
    conv_rad = np.zeros(chunks)

    # convolution mask array
    new_conv = np.zeros([chunks,no_cells,no_cells])
    
    # coarse mask
    scale_down = 2  # making the grid coarser    
    
    # store all profile masks
    fine_mask_len = profile.get_fine_mask_len(resolution)
    nbig = fine_mask_len*2
    nsmall = int(nbig/scale_down)
    addition_masks =np.zeros([chunks,nsmall,nsmall])
    
    # loops through the list of dataframes each ordered by ascending mass
    
    for j in range(0,chunks):
        Mvir_avg[j] = np.mean((df['Mvir'][bin_markers[j]:bin_markers[j+1]])) / h
        conv_rad[j] = (1+redshift)*((Mvir_avg[j])**(1/3) / Rvir_den(redshift)) # comoving radius
              
        # NFW 
        # add redshift to the function above    
        R_s= conv_rad[j]/(halo_conc(redshift,Mvir_avg[j])*cellsize)  
        rho_nought = rho_0(redshift,Mvir_avg[j],R_s)

        fine_mask = profile.get_mask(Mvir_avg[j], conv_rad[j], redshift, resolution, scaling_radius, cellsize)
        

        # spherical tophat
        if profile == 'tophat_spherical':
            r = (x**2+y**2)**.5
            fine_mask = r <= (scaling_radius*scale_down*conv_rad[j]/cellsize)
            
            fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
            Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
            fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
        elif profile == 'NFW':
            
            vec_integral = np.vectorize(NFW2D)
            fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)
            
            r=(x**2+y**2)**.5 # * scale_down
            
            fine_mask=fine_mask.astype(float)
            fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
            
       # elif profile == 'custom':
            
            # TODO the 'custom' profile I got from Adnan was broken when handed to me. Not sure what it was supposed to do. 

            # Functions for profile
            # Currently hardcoded but should allow the user to provide the function as inputs
            
            #f1_1 = lambda x,y,z: 1/(x**2+y**2+z**2+.5)**.5
            #f1_2 = lambda x,y,z: 1/(x**2+y**2+z**2+.5)
            
            # Radius of first profile:
            #R= (scaling_radius*scale_down*conv_rad[j]/cellsize)/2 
            
            #vec_integral = np.vectorize(func3Dto2D)
            
            #mask1 = x**2+y**2 < (R/2)**2
            #mask1 = mask1.astype(float)*vec_integral(f1_1,x,y,R/2)

            #mask2_1 = x**2+y**2 >= (R/2)**2
            #mask2_1= mask2_1.astype(float)
            #mask2_2 = x**2+y**2 <= (R)**2
            #mask2_2= mask2_2.astype(float)
            #mask2=mask2_1*mask2_2
            #mask2 = mask2*vec_integral(f1_2,x,y,R)

            #fine_mask=mask1+mask2
            
        #elif profile == 'custom_tophat':
            
            # TODO the 'custom_tophat' profile I got from Adnan was broken when handed to me. Not sure what it was supposed to do. 

            # Functions for profile
            # Currently hardcoded but should allow the user to provide the function as inputs
            
            # Radius of first profile:
            #R= (scaling_radius*scale_down*conv_rad[j]/cellsize)
            
            #f1_1 = lambda x,y,z: np.exp(-((x**2+y**2+z**2)/R**2)**30)
            
            #vec_integral = np.vectorize(func3Dto2D)
            
            #mask1 = vec_integral(f1_1,x,y,R)

            #fine_mask=mask1
        
        # testing code to add Fire simulation halos
        elif profile == "fire":
            Msun =  1.9889e33  # gr 

  
            # This is FIRE simulation profile (from eyeballing https://arxiv.org/pdf/1811.11753.pdf and using their $r^{-2}$ scaling)
            #  I'm using the $r^{-2}$ profile they find, with an exponetial cutoff at rmax, where we use conservation of mass to determine rmax.  So
            #
            #$$\rho = \rho_{0}\left(\frac{r}{r_*} \right)^{-2} \exp[-r/r_{\rm max}]$$
            #
            #They results technically hold for $10^{10} - 10^{12} M_\odot$ and $z=0$, but I think it's reasonable to assume they extrapolate to other moderate redshifts and somewhat more massive halos.
            #
            #To normalize things we are using that 
            #$$M_{\rm gas} = \int 4 \pi r^2 dr \rho_{0} \left(\frac{r}{r_*} \right)^{-2} \exp[-r/r_{\rm max}] = 4 \pi r_{\rm max} \rho_0 r_*^2$$
            #
            #Note that sometimes I'm working with number densities and sometimes mass densities
            #

            adjustmentfactor = 1  #probably we want to range between 0.5-2 as the numbers look sensible in this range
            
            #These are specifications taken from the fire simulations
            RinterpinRvir = 0.3  # this is the point where I read off the density nrmalization
            logMinterp = np.array([10., 11., 12.])  # these are log10 of the halo masses they consider
            nHinterp = np.array([0.5e-4, 0.8e-4, 1e-4])  #these are their number densities in cubic cm
                
            nHinterp = interp1d(logMinterp, nHinterp, fill_value="extrapolate")
            Mhalo = Mvir_avg[j]
            Rvir = conv_rad[j]

            rho0 = adjustmentfactor*nHinterp(np.log10(Mhalo))
            Rinterp = RinterpinRvir* Rvir #rvir(Mhalo, z)
            
            Ntot = Mhalo*Msun*fb/(mu*mp)/(MPCTOCM**3) # TODO msun here is is grams, but elsewhere it is in kg. Double check math.
            rmax = Ntot/(4.*np.pi*rho0*Rinterp**2 )  #from integrating above expression for rho
            
            # TODO these two lines are dead code, they were returned but not used before.
            rarr = np.logspace(-2, np.log10(5*rmax), 100) # Cut off at 5 times exponential cutoff
            rhoarr = rho0*(rarr/Rinterp)**-2*np.exp(-rarr/rmax) #number density: per cm3
            
            ## creating a mask
            # TODO perf bottleneck is here
            f1 = lambda x, y, z: fire_func(((x**2+y**2+z**2)**.5), rmax, Rinterp, rho0, cellsize)
            
            vec_integral = np.vectorize(fire3Dto2D)   

            mask1 = vec_integral(f1,x,y,rmax/cellsize)
            fine_mask = mask1.astype(float)
        
        # Precipitation model
        elif profile == "precipitation":       
            logMhalo = np.log10(Mvir_avg[j])
            Rvirkpc = 1000*conv_rad[j]
            XRvir = 2 # XRvir is how many virial radii to go out
            #Zmetal is the metalicity; tratcrit is the coolin crieteria -- both of these don't need to change
            Zmetal = 0.3
            tratcrit = 10
    
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
            reducedarr = fitarray[(fitarray[:, 3] == Zmetal) & (fitarray[:, 2] ==  tratcrit)]
            reducedarr = reducedarr[::-1] #reverses array
            
            #old interpolation
            #n1 = np.interp(logMhalo, np.log10(reducedarr[:, 1]), reducedarr[:, 6])
            #xi1 = np.interp(logMhalo, np.log10(reducedarr[:, 1]), reducedarr[:, 7])
            #n2 = np.interp(logMhalo, np.log10(reducedarr[:, 1]), reducedarr[:, 8])
            #xi2 = np.interp(logMhalo, np.log10(reducedarr[:, 1]), reducedarr[:, 9])                  
            
            #better interpolation
            logn1 = interp1d(np.log10(reducedarr[:, 1]), np.log10(reducedarr[:, 6]), kind='linear', fill_value='extrapolate')(logMhalo)   
            n1 = 10**logn1
            xi1 = interp1d(np.log10(reducedarr[:, 1]), reducedarr[:, 7], kind='linear', fill_value='extrapolate')(logMhalo) 
            logn2 = interp1d(np.log10(reducedarr[:, 1]), np.log10(reducedarr[:, 8]), kind='linear', fill_value='extrapolate')(logMhalo)   
            n2 = 10**logn2
            xi2 = interp1d(np.log10(reducedarr[:, 1]), reducedarr[:, 9], kind='linear', fill_value='extrapolate')(logMhalo) 
                
            #print(logMhalo, n1, xi1, n2, xi2)

            rkpc = np.logspace(0, np.log10(Rvirkpc*XRvir), 300)#number of radial bins

            #Voit 2018 fitting formulae
            rhoarr = np.array([rkpc, 1/np.sqrt(1/(n1*(rkpc)**-xi1+ 1e-20)**2 + 1/(n2*(rkpc/100)**-xi2 + 1e-20)**2)])
            
            #Integrate to see how much mass is missed by this profile  (I've checked these seems reasonable)
            rhointerp = interp1d(np.log(rhoarr[0]), 4.*np.pi*rhoarr[0]**3*rhoarr[1], kind='cubic', fill_value='extrapolate')
            conv = (1.989e33/(1.67e-24))/(3.086e21)**3  #constants: use your values TODO
            mtotal = integrate.quad(rhointerp, 0, 3*np.log(10))[0]/conv
            
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
            y,x = np.ogrid[-20*resolution: 20*resolution, -20*resolution:20*resolution]

        #     f1= lambda x, y, z: my_func(((x**2+y**2+z**2)**.5), n1,n2,xi1,xi2,neconstant,cellsize_kpc)
            f1= lambda x, y, z: precipitation_func(((x**2+y**2+z**2)**.5), n1,n2,xi1,xi2,neconstant,cellsize_kpc,Rvirkpc,XRvir,redshift)
                     
            vec_integral=np.vectorize(fire3Dto2D)
            
        #     mask1 = vec_integral(f1,x,y,XRvir*Rvirkpc/cellsize_kpc)
            mask1 = vec_integral(f1,x,y,XRvir*Rvirkpc/cellsize_kpc)
            r=(x**2+y**2)**.5 # * scale_down
            mask1=mask1.astype(float)
            mask1[r > (XRvir*Rvirkpc/cellsize_kpc)]=0         

        #     mask1[r <= (XRvir*Rvirkpc/cellsize_kpc)] =+ neconstant 
               
            fine_mask = mask1


        
        elif profile == "2RVSTH_and_NFW_13.5":
            if Mvir_avg[j] <= 10**13.5:
                
                r = (x**2+y**2)**.5
                
                # Don't hardcode scaling radius, change later
                fine_mask = r <= ((2*scaling_radius)*scale_down*conv_rad[j]/cellsize)
            
                fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
                Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
                fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
            elif Mvir_avg[j] > 10**13.5:
                vec_integral = np.vectorize(NFW2D)
                fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)

                r=(x**2+y**2)**.5 # * scale_down

                fine_mask=fine_mask.astype(float)
                fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
        
        elif profile == "2RVSTH_and_NFW_13":
            if Mvir_avg[j] <= 10**13:
                
                r = (x**2+y**2)**.5
                
                # Don't hardcode scaling radius, change later
                fine_mask = r <= ((2*scaling_radius)*scale_down*conv_rad[j]/cellsize)
            
                fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
                Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
                fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
            elif Mvir_avg[j] > 10**13:
                vec_integral = np.vectorize(NFW2D)
                fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)

                r=(x**2+y**2)**.5 # * scale_down

                fine_mask=fine_mask.astype(float)
                fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
        
        elif profile == "2RVSTH_and_NFW_12.5":
            if Mvir_avg[j] <= 10**12.5:
                
                r = (x**2+y**2)**.5
                
                # Don't hardcode scaling radius, change later
                fine_mask = r <= ((2*scaling_radius)*scale_down*conv_rad[j]/cellsize)
            
                fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
                Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
                fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
            elif Mvir_avg[j] > 10**12.5:
                vec_integral = np.vectorize(NFW2D)
                fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)

                r=(x**2+y**2)**.5 # * scale_down

                fine_mask=fine_mask.astype(float)
                fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
        
        
        elif profile == "2RVSTH_and_NFW_12":
            if Mvir_avg[j] <= 10**12:
                
                r = (x**2+y**2)**.5
                
                # Don't hardcode scaling radius, change later
                fine_mask = r <= ((2*scaling_radius)*scale_down*conv_rad[j]/cellsize)
            
                fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
                Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
                fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
            elif Mvir_avg[j] > 10**12:
                vec_integral = np.vectorize(NFW2D)
                fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)

                r=(x**2+y**2)**.5 # * scale_down

                fine_mask=fine_mask.astype(float)
                fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
        
        elif profile == "2RVSTH_and_NFW_11.5":
            if Mvir_avg[j] <= 10**11.5:
                
                r = (x**2+y**2)**.5
                
                # Don't hardcode scaling radius, change later
                fine_mask = r <= ((2*scaling_radius)*scale_down*conv_rad[j]/cellsize)
            
                fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
                Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
                fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
            elif Mvir_avg[j] > 10**11.5:
                vec_integral = np.vectorize(NFW2D)
                fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)

                r=(x**2+y**2)**.5 # * scale_down

                fine_mask=fine_mask.astype(float)
                fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
        
        elif profile == "2RVSTH_and_NFW_11":
            if Mvir_avg[j] <= 10**11:
                
                r = (x**2+y**2)**.5
                
                # Don't hardcode scaling radius, change later
                fine_mask = r <= ((2*scaling_radius)*scale_down*conv_rad[j]/cellsize)
            
                fine_mask=fine_mask.astype(float)
#             mask5 = mask5* (2*(((scaling_radius*scale_down*conv_rad[j]/cellsize)**2-(r**2))**2)**.25)
            
                Rv= (scaling_radius*scale_down*conv_rad[j]/cellsize)
                fine_mask = fine_mask* ((((1)**2-((r/Rv)**2))**2)**.25)
            
            elif Mvir_avg[j] > 10**11:
                vec_integral = np.vectorize(NFW2D)
                fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)

                r=(x**2+y**2)**.5 # * scale_down

                fine_mask=fine_mask.astype(float)
                fine_mask[r> scale_down*conv_rad[j]/cellsize] =0      
        else:
            raise ValueError("Not valid profile provided")
                
        
        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        nsmall = int(nbig/scale_down)
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)

        # Area of cells needed for normalization
        totalcellArea4 = 0
        totalcellArea4 = sum(sum(coarse_mask))* ((cellsize)**2)
        
        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])    
        
        # The coordinates are being multiplied by 4 to yield the halo coordinates on the 1024 grid
        ix = ((((np.around(4*resolution*((df[bin_markers[j]:bin_markers[j+1]]['x'].values)/(250/256))))))%(resolution*1024)).astype(int)
        iy = ((((np.around(4*resolution*((df[bin_markers[j]:bin_markers[j+1]]['y'].values)/(250/256))))))%(resolution*1024)).astype(int)  
        xy=(ix,iy)

        # issue: the method does not add repeated coordinates BUG is that right?
        halo_cell_pos[xy] += 1

        # convolve the mask and the halo positions
        new_conv[j,:,:] = (Mvir_avg[j]/(totalcellArea4)) * my_convolve(halo_cell_pos,coarse_mask)

        # store addition masks
        addition_masks[j,:,:]= (Mvir_avg[j]/(totalcellArea4))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM)*coarse_mask
        
    
    return (new_conv.sum(0))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM), conv_rad, addition_masks, Mvir_avg


# Halos removed field

#This function combines many steps of subtraction and smoothing to yield a density field from which halos have been removed
def halos_removed_field(current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,subtraction_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc):
    
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    bin_markers= halo_array_for_convolution[1]
    
    # convolve halos
    subtraction_profile = subtract_halos(df,resolution,bin_markers,subtraction_halo_profile,scaling_radius,redshift)[0]
    subtraction_profile_smooth = gauss_sinc_smoothing(subtraction_profile,sigma_gauss,width_sinc,1)
    
    # create coarse grid
    subtracted_coarse= smoothfield(subtraction_profile_smooth,1024,den_grid_size)
    
    # remove halos from the density field
    halos_removed_coarse=removeConvolvedHalos(density_field,subtracted_coarse)
    
    return halos_removed_coarse,subtracted_coarse


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
    sigma_gauss = 4*resolution/np.sqrt(3)

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
                           addition_profile,scaling_radius,resolution,sigma_gauss,width_sinc)
        
        
        halos_reAdded[i,:,:] = conv_all_steps[0]
        
        halos_subtraction_coarse[i,:,:] = halos_removed[1]
        halo_field[i,:,:] = conv_all_steps[2]
        halo_masks[i,:,:,:]= conv_all_steps[3]
        halos_removed_fields[i,:,:] =  halos_removed[0]
        
        
    # returns halo array and masks used to add halos back
    return halos_reAdded,halo_masks,halos_subtraction_coarse,halo_field,halos_removed_fields, conv_all_steps[5],conv_all_steps[6]



def hist_profile(sim_provider: SimulationProvider, den_grid_size, RS_array, min_mass, max_mass,
                                       log_bins, subtraction_halo_profile, addition_profile: CGMProfile, scaling_radius, resolution):
    """
    This function runs everything needed apply cgmbrush.

    Outputs: histograms, halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, 
    halos removed field, stacked halo field, virial radii, halo masses
    """
    
    # halo array
    t = halo_subtraction_addition(sim_provider,den_grid_size,RS_array,min_mass,max_mass,
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
    
    if len(RS_array)==1:
        t6 = t1
    
    else:
        t6 = stack_all_arrays(t1,RS_array)
        
    t7 = create_histograms(t6,resolution) # TODO turn off?
    
    t8 = t[5]
    t9 = t[6]
    
    # Outputs: 
    # histograms, halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, halos removed field, stacked halo field, virial radii, halo masses
    return t7,t1,t2,t3,t4,t5,t6,t8,t9



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

# Calculate Dispersion measure for each mass bin
def DM_for_mass_bin(halo_array,df,crop_grid):
    
    res=halo_array.shape[0]
    trim_dim=crop_grid
#     DM_rad = np.zeros([len(df), 2*trim_dim+1,2*trim_dim+1])
    
    DM_rad = np.zeros([len(df), trim_dim,trim_dim])
    
#     ix = (df['ix'].values)
#     iy = (df['iy'].values)
    
    ix = ((np.around((4*res/1024)*((df['x'].values)/(250/256))))%(res)).astype(int)
    iy = ((np.around((4*res/1024)*((df['y'].values)/(250/256))))%(res)).astype(int)
#     print(ix)
    
    for i in range(0,len(df)):            
            
        if ix[i] > int(trim_dim) and ix[i] < res - int(trim_dim) and iy[i] > int(trim_dim) and iy[i] < res- int(trim_dim):
#             trimmed = halo_array[ix[i]-int(trim_dim):ix[i]+int(trim_dim+1),iy[i]-int(trim_dim):iy[i]+int(trim_dim +1)]
            trimmed = halo_array[ix[i]-int(trim_dim/2):ix[i]+int(trim_dim/2),iy[i]-int(trim_dim/2):iy[i]+int(trim_dim/2)]
            DM_rad[i,:,:] = trimmed

    return DM_rad
        
# Function: calculates the radial profile for each mass bin
def radial_profile_array(DM_halo_array,halos_per_mass_bin,len_rad_ar):
    
    trim_dim=DM_halo_array.shape[1]
    center=int(trim_dim/2)
    
    DM_mass_bin= np.zeros([len(halos_per_mass_bin),trim_dim,trim_dim])
    rad_prof_mass_bin = np.zeros([len(halos_per_mass_bin),len_rad_ar])

    
    # The -1 has to be thought more about
    for i in range(0,len(halos_per_mass_bin)-1):
        DM_mass_bin[i,:,:]= np.mean(DM_halo_array[halos_per_mass_bin[i]:halos_per_mass_bin[i+1],:,:],axis=0)
        rad_prof_mass_bin[i,:] = radial_profile(DM_mass_bin[i,:,:],center,center)
    
    return rad_prof_mass_bin,DM_mass_bin


# Single function for radial profile of DM for a given DM array and grid size
def DM_vs_radius(DM_array,halo_data_frame,crop_grid_dim, mass_bins, hist_mass_bins):
    
#     len_rad_ar = radial_profile(DM_for_mass_bin(DM_array,halo_data_frame,crop_grid_dim)[0,:,:],int(crop_grid_dim),int(crop_grid_dim)).shape[0]
    
    # I can run profile of masks function to get this instead of running the above function
    len_rad_ar=15
    
    trimmed_ar = DM_for_mass_bin(DM_array,halo_data_frame,crop_grid_dim)
#     print((trimmed_ar.shape)[-1])
    
    final_ar=radial_profile(trimmed_ar[0,:,:],(trimmed_ar.shape)[-1]/2,(trimmed_ar.shape)[-1]/2)
        
#     return radial_profile_array(DM_for_mass_bin(DM_array,halo_data_frame,crop_grid_dim),hist_mass_bins,((final_ar.shape)[0]))

    output = radial_profile_array(trimmed_ar,hist_mass_bins,((final_ar.shape)[0]))
    
    return output[0],output[1],trimmed_ar
    
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
# def translate_array(array,seed_x,seed_y):
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

# Stack all arrays
# Translate and redshift arrays 
def stack_all_arrays(halos_reAdded,RS_array):
    
    resolution= np.shape(halos_reAdded)[-1]
    halos_reAdded_translated = np.zeros([len(RS_array),resolution,resolution])
    halos_reAdded_translated[0,:,:] = redshifted_DM(halos_reAdded[0,:,:],RS_array[0])
    
    for i in range(1, len(RS_array)):
        halos_reAdded_translated[i,:,:] = redshifted_DM(translate_array(halos_reAdded[i,:,:]),RS_array[i]) # 10,10 are seeds to generate random numbers
    
    
    return halos_reAdded_translated

def create_histograms(halos_reAdded_translated,resolution):
    nbin=80
    hist = histArray(sum(halos_reAdded_translated[:,:,:]),nbin,int(resolution),0,3*np.mean(sum(halos_reAdded_translated[:,:,:])))

    return hist


# Create histogram for the stack without min and max defined

# def create_histograms(halos_reAdded_translated,resolution):
#     nbin=1000
#     hist = histArray(sum(halos_reAdded_translated[:,:,:]),nbin,int(1024*resolution))

#     return hist






























###########################################
# Utils. These are to help run the code. 
# TODO Perhaps they will go in a different part of the package.
###########################################

varFolder = "../var"
testFolder = "test"

# Intermediate numpy arrays get can be saved into var folder outside version control
def saveFig(filename_base, fig):
    file_path = os.path.join(varFolder, filename_base)
    
    if not(os.path.exists(varFolder)):
        os.makedirs(varFolder)

    fig.savefig(file_path + '_images')


# Intermediate numpy arrays get can be saved into var folder outside version control
def saveArray(filename, *arrays, folder = varFolder):
    file_path = os.path.join(folder, filename)
    
    if not(os.path.exists(folder)):
        os.makedirs(folder)

    if len(arrays) == 1: # unwrap it out of tuple or it will not round-trip
        np.save(file_path, arrays[0]) 
    else:
        np.save(file_path, arrays) 
    
def loadArray(filename, folder = varFolder):
    file_path = os.path.join(folder, filename + ".npy")
    try:
        return np.load(file_path, allow_pickle=True)
    except FileNotFoundError:
        file_path = os.path.join(folder, filename + ".txt")
        return np.load(file_path, allow_pickle=True)


class Configuration:
    """Configuration for a run of cgmbrush."""

    # Default options
    def __init__(self, addition_halo_profile, scaling_radius, resolution=1, file_prefix=None, den_grid_size=256, RS_array=[0], load_from_files=False):
        
        # Profile to use for adding in CGM
        self.addition_halo_profile = addition_halo_profile
        self.file_prefix = file_prefix
        if (self.file_prefix == None):
            self.file_prefix = self.addition_halo_profile

        self.scaling_radius = scaling_radius

        # Resolution: choose between 256 and 512 grid
        if den_grid_size != 256 and den_grid_size != 512:
            raise ValueError("Only resolutions 256 and 512 are allowed")
        self.den_grid_size = den_grid_size 

        # User provides a redshift array
        self.RS_array = RS_array # For a single box, we only use the redshift 0 box

        # Profile used for subtracting halos from the density field
        self.subtraction_halo_profile = 'NFW'

        self.min_mass = 10**10 # TODO into config
        self.max_mass = 10**14.5
        self.log_bins = 30

        self.resolution = resolution # x1024

        self.load_from_files = load_from_files
        self.results = None
        self.figure = None
    

    def run(self, plots=False, trace=False):
        """Run this configuration."""

        filename = self.file_prefix + str(self.resolution) + '_' + str(self.den_grid_size) + "_" + str(datetime.date.today())

        if self.load_from_files:
            try:
                self.results = loadArray(filename)
            
            except IOError:
                print("Cache miss: " + filename)
                #pass # File cache doesn't exist, swallow and compute it instead

        if self.results == None:
            print("Performing Calculations for " + filename)
                                   
            if trace:
                pr = cProfile.Profile()
                pr.enable()

            self.results = hist_profile(self.provider, self.den_grid_size, self.RS_array, self.min_mass, 
                                                self.max_mass, self.log_bins, self.subtraction_halo_profile, 
                                                self.addition_halo_profile, self.scaling_radius, self.resolution)
            saveArray(filename, *self.results)
            
            if trace:
                pr.disable()
                s = io.StringIO()
                ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
                ps.print_stats()

                version = 1
                # TODO auto version incrementing
                perf_file_path = os.path.join(varFolder, 'perf_convo_' + filename + '_v%s.txt' % version)
                with open(perf_file_path, 'w') as f:
                    f.write(s.getvalue())

        if plots:
            original = self.provider.get_density_field(0, self.den_grid_size)
            background_dm = self.results[5][0]
            cgm_only = self.results[4][0]
            density_final = self.results[1][0]

            fig, axes = plt.subplots(2,2,figsize=(24, 24))
            plt.tight_layout()
            pos = axes[0][0].imshow(original) 
            fig.colorbar(pos, ax=axes[0][0])
            axes[0][0].title.set_text('Original Density Field')
            pos = axes[0][1].imshow(background_dm) 
            fig.colorbar(pos, ax=axes[0][1])
            axes[0][1].title.set_text('Density minus halos')
            pos = axes[1][0].imshow(cgm_only) 
            fig.colorbar(pos, ax=axes[1][0])
            axes[1][0].title.set_text('CGM Profile to add')
            pos = axes[1][1].imshow(density_final) 
            fig.colorbar(pos, ax=axes[1][1])
            axes[1][1].title.set_text('Final Product')

            self.figure = fig
            saveFig(filename, fig)