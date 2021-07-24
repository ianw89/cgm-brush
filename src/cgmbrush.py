from __future__ import print_function 
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from numpy import genfromtxt
from random import randrange
from numpy.polynomial import Polynomial as P
from IPython.display import Latex
import pandas as pd
import os

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