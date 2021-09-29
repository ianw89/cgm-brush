from __future__ import print_function 
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline
from matplotlib.patches import Rectangle
from cgmbrush import *

# Resolution: choose between 256 and 512 grid
den_grid_size = 256
# User provides a redshift array
RS_array = [0] # For a single box, we only use the redshift 0 box
# Profile used for subtracting halos from the density field
subtraction_halo_profile = 'NFW'
# Mass range of halos
min_mass=10**10
max_mass=10**14.5
log_bins=30

provider = BolshoiProvider()
orig_den_256 = provider.get_density_field(0, 256)
halos_0 = provider.get_halos(0)


STH8_256 = loadArray('STH8_256_2021-09-13') # for compatability with how this file was written

# This function orders the halo array in ascending mass
# Outputs: ordered dataframe of halos, mass bins
df= create_halo_array_for_convolution(halos_0,min_mass,max_mass,log_bins)

# Mean DM of single box
mean_DM=np.mean(orig_den_256)

# Mass bins out of the 30 bins
M_chosen = [1,10,12,18,25]

# Table of virial radii and avg masses. TODO do these vary from run to run?
vir_rad_ar = STH8_256[7]
avg_mass_ar = STH8_256[8]

# Specify resolution
resolution=8
grid_size = resolution*1024

# dimension of the small grid around the halo we want to crop
trim_dim=int((10*resolution))

# Radial extent of the plots in Mpc
extent = (L/grid_size)*(trim_dim/2)

folder = '/Volumes/Seagate Backup Plus Drive/CGM-FRB-Data/'
load_data = True
load_DM_vs_rad = True
load_masks = False

STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=resolution, folder=folder)
STH_config.datestamp = '2021-09-28'
STH_config.run(load_from_files=load_data)
STH_config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
STH_config.generate_profile_of_masks(load_from_files=load_masks)
STH_config.clear_results()
#STH_256 = STH_config.results_as_tuple # for compatability with how this file was written
STH8_DMvsR = STH_config.DM_vs_R1
STH8_masks = STH_config.mask_profiles


"""

STH_2_config = Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=resolution, folder=folder)
STH_2_config.run(load_from_files=True)
STH_2_config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
STH_2_config.generate_profile_of_masks(load_from_files=load_masks)
STH_2_config.clear_results()
#STH_2_256 = STH_2_config.results_as_tuple # for compatability with how this file was written

FIRE_config = Configuration(FireProfile(), 1, provider=provider, resolution=resolution, folder=folder)
FIRE_config.run(load_from_files=True)
FIRE_config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
FIRE_config.generate_profile_of_masks(load_from_files=load_masks)
FIRE_config.clear_results()
#FIRE_256 = FIRE_config.results_as_tuple # for compatability with how this file was written

NFW_config = Configuration(NFWProfile(), 1, provider=provider, resolution=resolution, folder=folder)
NFW_config.run(load_from_files=True)
NFW_config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
NFW_config.generate_profile_of_masks(load_from_files=load_masks)
NFW_config.clear_results()
#NFW_256 = NFW_config.results_as_tuple # for compatability with how this file was written

P_config = Configuration(PrecipitationProfile(), 1, provider=provider, resolution=resolution, folder=folder)
P_config.run(load_from_files=True)
P_config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
P_config.generate_profile_of_masks(load_from_files=load_masks)
P_config.clear_results()
#P_256 = P_config.results_as_tuple # for compatability with how this file was written

STH8_2_DMvsR = STH_2_config.DM_vs_R1
NFW8_DMvsR = NFW_config.DM_vs_R1
fire8_DMvsR = FIRE_config.DM_vs_R1
P8_DMvsR = P_config.DM_vs_R1

STH8_2_masks = STH_2_config.mask_profiles
NFW8_masks = STH_2_config.mask_profiles
fire8_masks = FIRE_config.mask_profiles
P8_masks = P_config.mask_profiles
"""
# BUG Kernal dies when generate_DM_vs_radius_profile() runs for res=8 (I think)


#####
# Plot generation
#####
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 30
XBIG_SIZE = 20
axis_fontsize = 20
curve_thickness = 2

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=XBIG_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


M_chosen = [1,10,15,20,25]

MpctoKpc =1000
mask_rescale = np.sqrt(2)
x_start=0  #kpc: min impact parameter 

DM_Rad_fig, DM_Rad_axs = plt.subplots(3, 1,
                        gridspec_kw={'hspace': 0.015, 'wspace': .2},figsize=(20,30))


# X-axis: mask grid has a diagonal of length sqrt(2) that needs to be factored in after the profile is calculated
x_axis = np.sqrt(2)*MpctoKpc*np.linspace(0,extent,STH8_DMvsR.shape[1])



# First plot

# DM_Rad_axs[0].semilogx(np.linspace(0,extent,TH_DMvsR.shape[1]),TH_DMvsR[M_chosen[0],:]-mean_DM,'-', label= '2D tophat')
DM_Rad_axs[0].semilogx(x_axis,fire8_DMvsR[M_chosen[1],:]-mean_DM,'-', label= 'FIRE',lw=5,color='green')
DM_Rad_axs[0].semilogx(x_axis,P8_DMvsR[M_chosen[1],:]-mean_DM,'-' ,label= 'Precipitation',lw=5,color='c')
DM_Rad_axs[0].semilogx(x_axis,NFW8_DMvsR[M_chosen[1],:]-mean_DM,'-' ,label= 'NFW',lw=5,color='blue')
DM_Rad_axs[0].semilogx(x_axis,STH8_DMvsR[M_chosen[1],:]-mean_DM,'-', label= '3D Tophat',lw=5,color= 'red')
DM_Rad_axs[0].semilogx(x_axis,STH8_2_DMvsR[M_chosen[1],:]-mean_DM,'-', label= '3D Tophat 2$R_{vir}$',lw=5,color= 'orange')


# Masks
DM_Rad_axs[0].semilogx(x_axis,STH8_masks[M_chosen[1],:],'--', lw=2,color='red')
DM_Rad_axs[0].semilogx(x_axis,STH8_2_masks[M_chosen[1],:],'--', lw=2,color='orange')
DM_Rad_axs[0].semilogx(x_axis,fire8_masks[M_chosen[1],:],'--',lw=2,color='green')
DM_Rad_axs[0].semilogx(x_axis,NFW8_masks[M_chosen[1],:],'--', lw=2,color='blue')
DM_Rad_axs[0].semilogx(x_axis,P8_masks[M_chosen[1],:],'--', lw=2,color='c')

# axis title
# DM_Rad_axs[0].set_title('Mass = %.1E' % Decimal(df[2][M_chosen[1]]),fontsize=14)

# Second plot

DM_Rad_axs[1].semilogx(x_axis,fire8_DMvsR[M_chosen[2],:]-mean_DM,'-' ,label= 'fire',lw=5,color='green')
DM_Rad_axs[1].semilogx(x_axis,P8_DMvsR[M_chosen[2],:]-mean_DM,'-' ,label= 'Precipitation',lw=5,color='c')
DM_Rad_axs[1].semilogx(x_axis,NFW8_DMvsR[M_chosen[2],:]-mean_DM,'-', label= 'NFW',lw=5,color='blue')
DM_Rad_axs[1].semilogx(x_axis,STH8_DMvsR[M_chosen[2],:]-mean_DM,'-', label= '3D tophat',lw=5,color='red')
DM_Rad_axs[1].semilogx(x_axis,STH8_2_DMvsR[M_chosen[2],:]-mean_DM,'-', label= '3D tophat',lw=5,color='orange')

# Masks
DM_Rad_axs[1].semilogx(x_axis,STH8_masks[M_chosen[2],:],'--', lw=1.5,color='red')
DM_Rad_axs[1].semilogx(x_axis,STH8_2_masks[M_chosen[2],:],'--', lw=1.5,color='orange')
DM_Rad_axs[1].semilogx(x_axis,fire8_masks[M_chosen[2],:],'--',lw=1.5,color='green')
DM_Rad_axs[1].semilogx(x_axis,NFW8_masks[M_chosen[2],:],'--', lw=1.5,color='blue')
DM_Rad_axs[1].semilogx(x_axis,P8_masks[M_chosen[2],:],'--', lw=2,color='c')


# Error bar
# STH
# DM_Rad_axs[1].semilogx(r_star_ar*1000,error_1000_MW_STH,ls='--',drawstyle='steps',color='red',lw=.5)
# DM_Rad_axs[1].semilogx(r_star_ar*1000,error_1000_MW_STH2,ls='--',drawstyle='steps',color='orange',lw=.5)
# DM_Rad_axs[1].semilogx(r_star_ar*1000,error_1000_MW_fire,ls='--',drawstyle='steps',color='green',lw=.5)
# DM_Rad_axs[1].semilogx(r_star_ar*1000,error_1000_MW_NFW,ls='--',drawstyle='steps',color='blue',lw=.5)


# DM_Rad_axs[1].set_title('Mass = %.1E' % Decimal(df[2][M_chosen[2]]),fontsize=14)


# Third plot
DM_Rad_axs[2].semilogx(x_axis,fire8_DMvsR[M_chosen[3],:]-mean_DM,'-', label= 'Fire',lw=5,color='green')
DM_Rad_axs[2].semilogx(x_axis,P8_DMvsR[M_chosen[3],:]-mean_DM,'-' ,label= 'Precipitation',lw=5,color='c')
DM_Rad_axs[2].semilogx(x_axis,NFW8_DMvsR[M_chosen[3],:]-mean_DM, '-',label= 'NFW',lw=5,color='blue')
DM_Rad_axs[2].semilogx(x_axis,STH8_DMvsR[M_chosen[3],:]-mean_DM,'-', label= '3D tophat',lw=5,color='red')
DM_Rad_axs[2].semilogx(x_axis,STH8_2_DMvsR[M_chosen[3],:]-mean_DM,'-', label= '3D tophat',lw=5,color='orange')

DM_Rad_axs[2].semilogx(x_axis,STH8_masks[M_chosen[3],:],'--', lw=2,color='red')
DM_Rad_axs[2].semilogx(x_axis,STH8_2_masks[M_chosen[3],:],'--', lw=2,color='orange')
DM_Rad_axs[2].semilogx(x_axis,fire8_masks[M_chosen[3],:],'--',lw=2,color='green')
DM_Rad_axs[2].semilogx(x_axis,NFW8_masks[M_chosen[3],:], '--',lw=2,color='blue')
DM_Rad_axs[2].semilogx(x_axis,P8_masks[M_chosen[3],:],'--', lw=2,color='c')

# DM_Rad_axs[2].set_title('Mass = %.1E' % Decimal(df[2][M_chosen[3]]),fontsize=14)


# ticks

# ax.xaxis.grid(True, which='minor')
# DM_Rad_axs[0].xaxis.grid(axis='x', which='minor', bottom=True)

# DM_Rad_axs[0].tick_params(axis='x', which='minor', bottom=True)
# DM_Rad_axs[1].tick_params(axis='x', which='minor', bottom=True)

# DM_Rad_axs[0].tick_params(axis='x', which='minor',fontsize=10)
# DM_Rad_axs[2].xaxis.set_tick_params(width=5)



DM_Rad_axs[0].tick_params('both', length=10, width=4, which='major')
DM_Rad_axs[0].tick_params('both', length=10, width=4, which='minor')

DM_Rad_axs[1].tick_params('both', length=10, width=4, which='major')
DM_Rad_axs[1].tick_params('both', length=10, width=4, which='minor')


DM_Rad_axs[2].tick_params('both', length=10, width=4, which='major')
DM_Rad_axs[2].tick_params('both', length=10, width=4, which='minor')


# legend
DM_Rad_axs[0].legend(loc='right',prop={'size':28}, frameon=False)
# DM_Rad_axs[0].set_ylabel('DM - <DM> [pc $cm^{-3}$]',fontsize=30)
DM_Rad_axs[1].set_ylabel('DM - <DM> [pc cm$^{-3}$]',fontsize=50)
# DM_Rad_axs[2].set_ylabel('DM - <DM> [pc cm$^{-3}$]',fontsize=30)
DM_Rad_axs[2].set_xlabel('Impact Parameter [kpc]',fontsize=50)

# Adding virial radii
DM_Rad_axs[0].axvline(MpctoKpc*(vir_rad_ar[M_chosen[1]]), color='k', linestyle='--', linewidth=1)
DM_Rad_axs[1].axvline(MpctoKpc*(vir_rad_ar[M_chosen[2]]), color='k', linestyle='--', linewidth=1)
DM_Rad_axs[2].axvline(MpctoKpc*(vir_rad_ar[M_chosen[3]]), color='k', linestyle='--', linewidth=1)

# x and y axis limits
DM_Rad_axs[0].set_ylim([0, 130])
DM_Rad_axs[0].set_xlim([x_start, .7*MpctoKpc])

DM_Rad_axs[1].set_ylim(ymin=0,ymax=399)
DM_Rad_axs[1].set_xlim([x_start, .7*MpctoKpc])

DM_Rad_axs[2].set_ylim(ymin=0,ymax=1499)
DM_Rad_axs[2].set_xlim([x_start, .7*MpctoKpc])


# Rectangular patch 
DM_Rad_axs[0].add_patch(Rectangle((0,0), 45, 130,facecolor='yellow'))
DM_Rad_axs[1].add_patch(Rectangle((0,0), 45, 399,facecolor='yellow'))
DM_Rad_axs[2].add_patch(Rectangle((0,0), 45, 1500,facecolor='yellow'))

# mass labels
DM_Rad_axs[0].text(310, 120, r'$6 \times 10^{11} M_\odot$ ',fontsize=34) #'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[1]]),fontsize=30)
DM_Rad_axs[1].text(310, 360,r'$3 \times 10^{12} M_\odot$ ',fontsize=34) #'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[2]]),fontsize=30)
DM_Rad_axs[2].text(310, 1350,r'$2 \times 10^{13} M_\odot$ ',fontsize=34) #'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[3]]),fontsize=30)

# DM_Rad_axs[0].text(2*MpctoKpc*.1, 65, 'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[1]]),fontsize=30)
# DM_Rad_axs[1].text(2*MpctoKpc*.1, 180, 'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[2]]),fontsize=30)
# DM_Rad_axs[2].text(2*MpctoKpc*.1, 925, 'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[3]]),fontsize=30)

# DM_Rad_axs[0].rc('xtick', labelsize=35)    # fontsize of the tick labels
# # plt.rc('ytick', labelsize=35)    # fontsize of the tick labels



DM_Rad_fig.savefig('DMvsRad_profiles_%s.pdf' % resolution, bbox_inches='tight')