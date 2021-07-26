from __future__ import print_function 
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from numpy import genfromtxt
from random import randrange
import random
import time
from numpy.polynomial import Polynomial as P
from IPython.display import Latex
import pandas as pd
import os
from scipy.ndimage.filters import convolve

########################################
# Paramters, Constants, and Functions
########################################

# Constants
Mpc = 3.096*10**24 # cm
msun=1.9891 * 10**30 #kilograms
mprot= 1.6726219 * 10**-27 #kilograms
Yhe =.25
nPS=msun/mprot*(1-Yhe/2) # accounts for the fact that some mass is in helium
lightspeed = 3e5 #km/s
pi=np.pi

# Cosmological parameters
z=0
h=.7
rho_c = 9.31000324385361 *10**(-30) * (Mpc**3 / (msun*1000))*(h/.7)**2 # M./Mpc^3 from #pcrit = 9.31000324385361e-30 g / cm3
OmegaM = 0.27
OmegaL = 1- OmegaM #assumes flat universe
OmegaB= 0.0469 
rho_m=OmegaM*rho_c
fd = 1 # fraction of baryons in diffused ionized gas (FRB paper)

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

def rho_vir(z):
    return (18*pi**2 - 82*q(z) - 39*q(z)**2)*(rho_c*(OmegaL + OmegaM *(1+z)**3))


def Rvir_den(z):
    return (4/3 * np.pi * rho_vir(z))**(1/3) # physical, units in 1/r**3

# q = OmegaL/ ((OmegaM*(1+z)**3)+ OmegaL)  
# rho_vir = (18*pi**2 - 82*q - 39*q**2)*(rho_c*(OmegaL + OmegaM *(1+z)**3))

#print(18*pi**2 - 82*q - 39*q**2, z)
# Rvir_den = (4/3 * np.pi * rho_vir)**(1/3) # physical, units in 1/r**3




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
def haloArray_minmax(pdHalosN,min, max):
    pdHalosN=pdHalosN.drop(pdHalosN[pdHalosN.Mvir < min].index)
    pdHalosN=pdHalosN.drop(pdHalosN[pdHalosN.Mvir > max].index)
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
# This section extracts the density field of the 256 grid from the simulation tables
########################################


# Function to import density and halo tables for a given redshift
# TODO I think this should be made to work with 512 resolution too
def import_density_field_256(redshift):
    
    # Have to hardcode z=0 table because of the unique column names
    if redshift == 0:
        # reading density field and halos data
        file_path= os.path.join(sims_folder, 'dens256-z-0.csv.gz')
        pdDens=pd.read_csv(file_path)

        # extracting columns
        pdDensN=pdDens[['Bolshoi__Dens256_z0__ix','Bolshoi__Dens256_z0__iy','Bolshoi__Dens256_z0__iz','Bolshoi__Dens256_z0__dens']]

        # 3D density array
        pdDensN=pdDensN.sort_values(['Bolshoi__Dens256_z0__ix','Bolshoi__Dens256_z0__iy','Bolshoi__Dens256_z0__iz'])
        tden = pdDensN['Bolshoi__Dens256_z0__dens'].values
        tden2=np.reshape(tden,(256,256,256))

        return ((tden2+1).sum(2))*10**6* dx* elecD(z) /(1+z)**2

    else:
        name = 'dens256-z-0.'+str(redshift)+'.csv.gz'
        den = pd.read_csv(name)
        den2=den[['Bolshoi__Dens256__ix','Bolshoi__Dens256__iy','Bolshoi__Dens256__iz','Bolshoi__Dens256__dens']]

        # 3D density array
        den_sorted=den2.sort_values(['Bolshoi__Dens256__ix','Bolshoi__Dens256__iy','Bolshoi__Dens256__iz'])
        den_vals = den_sorted['Bolshoi__Dens256__dens'].values
        den256=np.reshape(den_vals,(256,256,256))
        
        return normDM((den256+1).sum(2),0)
    #     return normDM((den256+1).sum(random.randint(0,1)),0)


# Create a single array of density fields for various redshifts
# Density fields: Choose 256 or 512 field
def extract_all_den_fields(RS_array,den_grid_size):

    all_den_fields = np.zeros([len(RS_array),den_grid_size,den_grid_size])

    if den_grid_size == 256:
        all_den_fields[0,:,:] = import_density_field_256(0)
    elif den_grid_size == 512:
        all_den_fields[0,:,:] = import_density_field_512(0) # BUG This is undefined... generalize the above?

    for i in range(1, len(RS_array)):
        if den_grid_size == 256:
            all_den_fields[i,:,:]=import_density_field_256(i)
        elif den_grid_size == 512:
            all_den_fields[i,:,:]=import_density_field_512(i)
    
    return all_den_fields 
    
    
# Extract halos for a given redshift
def extract_halos(redshift):

    name = 'halo-z-0.'+str(redshift)+'.csv.gz'    
    file_path= os.path.join(sims_folder, name)
    halos = pd.read_csv(file_path)
    return halos




########################################
# Convolution Functions
########################################


# Create halo array from halo table for convolution
def create_halo_array_for_convolution(pdHalos,M_min,M_max,logchunks):
    halos = haloArray_minmax(pdHalos,M_min,M_max)
    df = halos.sort_values(by='Mvir',ascending=True)
    sorted_haloMasses = df['Mvir'].values
    histH, binsH = np.histogram(sorted_haloMasses,bins=np.logspace(np.log10(M_min),np.log10(M_max), logchunks))
    binz=np.append(0,histH.cumsum())
    
    return df,binz,binsH

# Convert custom 3D function to 2D by integrating over z
def func3Dto2D(f,x,y,Rvir):
    return integrate.quad(f,-(Rvir**2-(x**2+y**2))**.5,(Rvir**2-(x**2+y**2))**.5,args=(x,y))[0]

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
# haloArray: array of halos
# resolution: 1024 or 2048
# chunks: no of chunks the halo array has been divided into -1
# bins: array cumulative sum of halos in different mass bins
# profile: tophat, NFW etc
# scaling_radius: scale radius for tophat halos
def subtract_halos(haloArray,mass_binz,resolution,chunks,bins,profile,scaling_radius,redshift):
    
    df = haloArray
    no_cells = 1024
    cellsize = L/(1024) 
    
    # array of halo masses and radii
    Mvir_avg = np.zeros(chunks)
    conv_rad = np.zeros(chunks)
    Rvir_avg = np.zeros(chunks)

    # convolution mask array
    convolution =np.zeros([chunks,no_cells,no_cells])
    
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
        Mvir_avg[j] = np.mean((df['Mvir'][bins[j]:bins[j+1]]))/h
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
        ix = ((((np.around(4*((df[bins[j]:bins[j+1]]['x'].values)/(250/256))))))%(1024)).astype(int)
        iy = ((((np.around(4*((df[bins[j]:bins[j+1]]['y'].values)/(250/256))))))%(1024)).astype(int)
        
        xy=(ix,iy)

        # issue: the method does not add repeated coordinates
        halo_cell_pos[xy] += 1
        
        # convolve the mask and the halo positions
        convolution[j,:,:] = (Mvir_avg[j]/(totalcellArea4))*convolve(halo_cell_pos,coarse_mask)
        
        
    
    return (convolution.sum(0))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM), conv_rad, Rvir_avg, fine_mask, coarse_mask,halo_cell_pos

    
# Status: This is the latest convolution function

# The function add halos

# arguments: 
# haloArray: array of halos
# resolution: 1024 or 2048
# chunks: no of chunks the halo array has been divided into -1
# bins: array cumulative sum of halos in different mass bins
# profile: tophat, NFW etc
# scaling_radius: scale radius for tophat halos

def add_halos(haloArray,mass_binz,resolution,chunks,bins,profile,scaling_radius,redshift):
    
    df = haloArray
    no_cells = 1024* resolution
    cellsize = L/(1024*resolution) 
    
    # array of halo masses and radii
    Mvir_avg = np.zeros(chunks)
    conv_rad = np.zeros(chunks)
    Rvir_avg = np.zeros(chunks)

    # convolution mask array
    convolution =np.zeros([chunks,no_cells,no_cells])
    
    # creating a coarse map out of a fine mask
    
    # fine mask
    # fine mask size has to correspond to the size of the mask that I eventually trim
    fine_mask_len = 20*resolution  
    fine_lower= -1*fine_mask_len
    fine_upper= fine_mask_len
    y,x = np.ogrid[fine_lower: fine_upper, fine_lower:fine_upper]
    
    # coarse mask
    scale_down = 2  # making the grid coarser
    smooth_r=2
    
    
#     coarse_mask_len = int(fine_mask_len/scale_down)
    fine_mask= np.zeros([fine_mask_len,fine_mask_len])
# #     fine_mask= np.zeros([2*coarse_mask_len,2*coarse_mask_len])
    
    
    # store all profile masks
    nbig = fine_mask_len*2
    nsmall = int(nbig/scale_down)
    addition_masks =np.zeros([chunks,nsmall,nsmall])
    
    # loops through the list of dataframes each ordered by ascending mass
    
    for j in range(0,chunks):
        Mvir_avg[j] = np.mean((df['Mvir'][bins[j]:bins[j+1]]))/h
        conv_rad[j] = (1+redshift)*((Mvir_avg[j])**(1/3) / Rvir_den(redshift)) # comoving radius
      
        # NFW 
        # add redshift to the function above
        R_s=0
        rho_nought=0
        
        R_s= conv_rad[j]/(halo_conc(redshift,Mvir_avg[j])*cellsize)  
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
            fine_mask =vec_integral(x,y,rho_nought,R_s,conv_rad[j]/cellsize)
            
            r=(x**2+y**2)**.5 # * scale_down
            
            fine_mask=fine_mask.astype(float)
            fine_mask[r> scale_down*conv_rad[j]/cellsize] =0 
            
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
        
        ## testing code to add Fire simulation halos
        elif profile == "fire":
        
            %run gasProfile_Fire.ipynb
            fine_mask = rhogasFire(Mvir_avg[j],conv_rad[j], redshift, adjustmentfactor,resolution)[2]
        
        # Precipitation model
        elif profile == "precipitation":
        
            %run gasProfile_precipitation.ipynb
            
#             fine_mask = rhogasFire(Mvir_avg[j],conv_rad[j], redshift, adjustmentfactor,resolution)[2]
            fine_mask = nePercipitation(np.log10(Mvir_avg[j]),1000*conv_rad[j],resolution,redshift)[1]
        
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
                
        
        # Smoothing method: reshaping
        # Generating coarse grid from fine grid: reshape method
        nsmall = int(nbig/scale_down)
        coarse_mask = fine_mask.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
        
        # Area of cells needed for normalization
        totalcellArea4=0
        totalcellArea4 = sum(sum(coarse_mask))* ((cellsize)**2)
        
        # populate array with halos
        halo_cell_pos = np.zeros([no_cells,no_cells])    
        
        # The coordinates are being multiplied by 4 to yield the halo coordinates on the 1024 grid
        ix = ((((np.around(4*resolution*((df[bins[j]:bins[j+1]]['x'].values)/(250/256))))))%(resolution*1024)).astype(int)
        iy = ((((np.around(4*resolution*((df[bins[j]:bins[j+1]]['y'].values)/(250/256))))))%(resolution*1024)).astype(int)
        
        
        xy=(ix,iy)

        # issue: the method does not add repeated coordinates
        halo_cell_pos[xy] += 1
        
        # convolve the mask and the halo positions
        convolution[j,:,:] = (Mvir_avg[j]/(totalcellArea4))*convolve(halo_cell_pos,coarse_mask)
        
        # store addition masks
        addition_masks[j,:,:]= (Mvir_avg[j]/(totalcellArea4))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM)*coarse_mask
        

    
    return (convolution.sum(0))*(Mpc**-3 *10**6)*nPS*(OmegaB/OmegaM), conv_rad, Rvir_avg, fine_mask, coarse_mask,halo_cell_pos,addition_masks,Mvir_avg


# Halos removed field

#This function combines many steps of subtraction and smoothing to yield a density field from which halos have been removed
def halos_removed_field(current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,subtraction_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc):
    
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    binz= halo_array_for_convolution[1]
    mass_binz = halo_array_for_convolution[2]
    
    # convolve halos
    subtraction_profile = subtract_halos(df,mass_binz,resolution,len(binz)-1,binz,subtraction_halo_profile,scaling_radius,redshift)[0]
    subtraction_profile_smooth = gauss_sinc_smoothing(subtraction_profile,sigma_gauss,width_sinc,1)
    
    # create coarse grid
    subtracted_coarse= smoothfield(subtraction_profile_smooth,1024,den_grid_size)
    
    # remove halos from the density field
    halos_removed_coarse=removeConvolvedHalos(density_field,subtracted_coarse)
    
    return halos_removed_coarse,subtracted_coarse


# Function subtracts and adds halos
def convolution_all_steps_final(current_halo_file,min_mass,max_mass,density_field,den_grid_size,redshift,log_bins,halos_removed_coarse,
                       addition_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc):
    t1 = time.time()
    
    # setup inputs for convolution
    
    halo_array_for_convolution = create_halo_array_for_convolution(current_halo_file,min_mass,max_mass,log_bins)
    df= halo_array_for_convolution[0]
    binz= halo_array_for_convolution[1]
    mass_binz = halo_array_for_convolution[2]
    
    t2 = time.time()
    
    # convolve halos for adding back
    addition_profile_initial=add_halos(df,mass_binz,resolution,len(binz)-1,binz,addition_halo_profile,scaling_radius,redshift)
    addition_profile = addition_profile_initial[0]
    addition_profile_masks=addition_profile_initial[6]
    
    t7 = time.time()
    
    # add halos to the subtracted field
    halosremoved_fine = (np.repeat((np.repeat(halos_removed_coarse,(1024/den_grid_size)*resolution,axis=0)),(1024/den_grid_size)*resolution,axis=1))
    roll_by = int((addition_profile.shape[0]/den_grid_size)/2)
    halosremoved_fine = np.roll(halosremoved_fine, -1*roll_by, axis=0)
    halosremoved_fine = np.roll(halosremoved_fine, -1*roll_by, axis=1)
    halos_added = addition_profile +  halosremoved_fine
    
    t8 = time.time()
    
    virial_rad = addition_profile_initial[1]
    halo_masses = addition_profile_initial[7]
    
    return halos_added, halosremoved_fine,addition_profile,addition_profile_masks,halos_removed_coarse,virial_rad,halo_masses
    


# Multiple Redshift Convolution
def halo_subtraction_addition(all_den_fields,den_grid_size,RS_array,min_mass,max_mass,log_bins,subtraction_halo_profile,
                             addition_halo_profile,scaling_radius,resolution):

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
#         current_halo_file = extract_halos(RS_array[i])
        current_halo_file = extract_halos(i)
        halos_removed = halos_removed_field(current_halo_file,min_mass,max_mass,all_den_fields[i,:,:],den_grid_size,RS_array[i],log_bins,subtraction_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc)
        
                
        conv_all_steps = convolution_all_steps_final(current_halo_file,min_mass,max_mass,all_den_fields[i,:,:],den_grid_size,RS_array[i],log_bins,halos_removed[0],
                           addition_halo_profile,scaling_radius,resolution,sigma_gauss,width_sinc)
        
        
        halos_reAdded[i,:,:] = conv_all_steps[0]
        
        halos_subtraction_coarse[i,:,:] = halos_removed[1]
        halo_field[i,:,:] = conv_all_steps[2]
        halo_masks[i,:,:,:]= conv_all_steps[3]
        halos_removed_fields[i,:,:] =  halos_removed[0]
        
        
    # returns halo array and masks used to add halos back
    return halos_reAdded,halo_masks,halos_subtraction_coarse,halo_field,halos_removed_fields, conv_all_steps[5],conv_all_steps[6]



# This function runs all the functions above functions at once
#
# Outputs:histograms,halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, 
         #halos removed field, stacked halo field
def hist_profile(all_den_fields,den_grid_size,RS_array,min_mass,max_mass,
                                       log_bins,subtraction_halo_profile,addition_halo_profile,scaling_radius,resolution):
    
    
    # halo array
    t = halo_subtraction_addition(all_den_fields,den_grid_size,RS_array,min_mass,max_mass,
                                       log_bins,subtraction_halo_profile,addition_halo_profile,scaling_radius,resolution)
    # Halos-readded field
    t1=t[0]
    
    # Halo addition masks
    t2=t[1]
    
    # Halos subtraction coarse
    t3=t[2]
    
    # Halo addition field
    t4=t[3]
    
    # Halos removed field
    t5 = t[4]
    
    if len(RS_array)==1:
        t6=t1
    
    else:
        t6= stack_all_arrays(t1,RS_array)
        
    t7= create_histograms(t6,resolution)
    
    t8 = t[5]
    t9 = t[6]
    
    # Outputs: 
    # 1s,halos-readded field, halo addition masks, halos subtraction coarse, halo addition field, halos removed field, stacked halo field
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
    
# BUG Mean DM of single box is calculated in MainFile notebook but used here in the functions below
#mean_DM=np.mean(import_density_field_256(0))
# TODO These 3 functions have no callers right now? Dead code?

# Function: plots radial profiles for different halo mass bins
def plot_radial_profile(radial_profile_array,mass_array,radial_extent,num_of_plots):
    
    plt.figure(figsize=(10,10))
    plot_x_dim = radial_profile_array.shape[1]
    
    for i in range(0,len(mass_array)-1,num_of_plots):
        
        plt.semilogx(np.linspace(0,radial_extent,plot_x_dim),radial_profile_array[i,:]-mean_DM, label= 'Mass = %.1E' % Decimal(mass_array[i]))
        plt.legend(loc='upper right')
        plt.xlabel('Radial distance (Mpc)')
        plt.ylabel('DM - <DM>')
    plt.savefig('DMvsRad_ImpactParam')

# Function: plots radial profiles for different halo mass bins
def plot_radial_profile_2(radial_profile_array_1,radial_profile_array_2,mass_array,radial_extent,num_of_plots):
    
    plt.figure(figsize=(10,10))
    plot_x_dim_1 = radial_profile_array_1.shape[1]
    print(plot_x_dim_1)
    plot_x_dim_2 = radial_profile_array_2.shape[1]
    print(plot_x_dim_2)
    
    for i in range(0,len(mass_array)-1,num_of_plots):
        
        plt.semilogx(np.linspace(0,radial_extent,plot_x_dim_1),radial_profile_array_1[i,:]-mean_DM, label= 'Mass = %.1E' % Decimal(mass_array[i]))
        plt.semilogx(np.linspace(0,radial_extent,plot_x_dim_2),radial_profile_array_2[i,:], linestyle='--',label= 'Mass = %.1E' % Decimal(mass_array[i]))
        plt.legend(loc='upper right')
        plt.xlabel('Radial distance (Mpc)')
        plt.ylabel('DM - <DM>')
        plt.ylim(0,250)
    plt.savefig('DMvsRad_ImpactParam')

def plot_radial_profile_3(radial_profile_array_1,radial_profile_array_2,radial_profile_array_3,mass_array,radial_extent,num_of_plots):
    
    plt.figure(figsize=(10,10))
    plot_x_dim_1 = radial_profile_array_1.shape[1]
    plot_x_dim_2 = radial_profile_array_2.shape[1]
    plot_x_dim_3 = radial_profile_array_3.shape[1]
    
    for i in range(0,len(mass_array)-1,num_of_plots):
        
        plt.semilogx(np.linspace(0,radial_extent,plot_x_dim_1),radial_profile_array_1[i,:]-75, label= 'Mass = %.1E' % Decimal(mass_array[i]))
        plt.semilogx(np.linspace(0,radial_extent,plot_x_dim_2),radial_profile_array_2[i,:]-75, linestyle='--',label= 'Mass = %.1E' % Decimal(mass_array[i]))
        plt.semilogx(np.linspace(0,radial_extent,plot_x_dim_3),radial_profile_array_3[i,:]-75, linestyle='-.',label= 'Mass = %.1E' % Decimal(mass_array[i]))
        plt.legend(loc='upper right')
        plt.xlabel('Radial distance (Mpc)')
        plt.ylabel('DM - <DM>')
        plt.ylim(0,300)
    plt.savefig('DMvsRad_ImpactParam')

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