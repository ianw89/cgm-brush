###################################################################################################
#
# DM_vs_R_plot.py 	                (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	        ianw89@live.com
#
###################################################################################################
from cgmbrush.cgmbrush import *
#halo = halo(cosmo)  #initialize halo classes  (probably more elegant way to do this...Matt: come back here)


#from cgmbrush import *
import numpy as np
import plotting_routines as makefig

M_chosen = [4,9,14]
provider = BolshoiProvider()

# Specify resolution
resolution =  8
grid_size = resolution*1024

load_data = False
load_DM_vs_rad = False
load_masks = False

series = []
date = '2021-3-14'
fire_date = '2021-3-14'


#print(cosmo.fb)
config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=1)
config.datestamp = date
vir_rad_ar = config.get_virial_radii()
avg_mass_ar = config.get_halo_masses()

print ("Masses chosen: ")
print ("{:e}".format(avg_mass_ar[M_chosen[0]]))
print ("{:e}".format(avg_mass_ar[M_chosen[1]]))
print ("{:e}".format(avg_mass_ar[M_chosen[2]]))

config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=resolution)
config.datestamp = date
config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
config.generate_profile_of_masks(load_from_files=load_masks)
STH_DMvsR = config.DM_vs_R1
STH_masks = config.mask_profiles
config.clear_results()
series.append((config.DM_vs_R1,config.mask_profiles,'3D Tophat','red'))

config = Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=resolution)
config.datestamp = date
config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
config.generate_profile_of_masks(load_from_files=load_masks)
STH_2_DMvsR = config.DM_vs_R1
STH_2_masks = config.mask_profiles
config.clear_results()
series.append((config.DM_vs_R1, config.mask_profiles, '3D Tophat 2$R_{vir}$', 'orange'))

config = Configuration(FireProfile(), 1, provider=provider, resolution=resolution)
config.datestamp = fire_date
config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
config.generate_profile_of_masks(load_from_files=load_masks)
fire_DMvsR = config.DM_vs_R1
fire_masks = config.mask_profiles
config.clear_results()
series.append((config.DM_vs_R1, config.mask_profiles, 'FIRE', 'green'))

config = Configuration(NFWProfile(), 1, provider=provider, resolution=resolution)
config.datestamp = date
config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
config.generate_profile_of_masks(load_from_files=load_masks)
NFW_DMvsR = config.DM_vs_R1
NFW_masks = config.mask_profiles
config.clear_results()
series.append((config.DM_vs_R1, config.mask_profiles, 'NFW', 'blue'))


config = Configuration(PrecipitationProfile(), 1, provider=provider, resolution=resolution)
config.datestamp = date
config.generate_DM_vs_radius_profile(load_from_files=load_DM_vs_rad)
config.generate_profile_of_masks(load_from_files=load_masks)
P_DMvsR = config.DM_vs_R1
P_masks = config.mask_profiles
config.clear_results()
series.append((config.DM_vs_R1, config.mask_profiles, 'Precipitation', 'c'))


redshift = 0
df = provider.get_halos(redshift)
df, bin_markers = create_halo_array_for_convolution(df, 10**10, 9*10**15, 30)



print("r200 = ", halo.r200Mz(cosmo, 1e12, 0))



fp = FireProfile()
    
for j in range(0,len(bin_markers) -1):
    Mvir_avg = np.mean((df['Mvir'][bin_markers[j]:bin_markers[j+1]])) / cosmo.h
    conv_rad = halo.comoving_rvir(cosmo, Mvir_avg, redshift) # comoving radius


    rvals, density = fp.get_analytic_profile(Mvir_avg, conv_rad)

    if j in M_chosen:
        plt.loglog(rvals, density)


# Radial extent of the plots in Mpc
# dimension of the small grid around the halo we want to crop
trim_dim=int((10*resolution))
dx = (provider.Lbox/grid_size)   
extent = (provider.Lbox/grid_size)*(trim_dim/2) 
MpctoKpc =1000
massbin = M_chosen[1]
#xmin = 0.5*extent/series[0][0].shape[1] #half a cell
#x_axis = np.sqrt(2)*MpctoKpc*np.linspace(xmin, extent-xmin, series[0][0].shape[1])
x_axis = MpctoKpc*(0.5+np.arange(0, series[0][0].shape[1]))*dx     #np.linspace(0.5*dx, extent-xmin, series[0][0].shape[1])

for data in series:
    Mass= np.sum(data[1][massbin,1:]*(x_axis[1:]**2-x_axis[0:-1]**2)*np.pi)
    Mass = Mass/cosmo.fb*Mpc/MpcInPc*mean_molecular_weight_electrons*mprot/msun*Mpc*Mpc/MpctoKpc**2
    print("label = ", data[2], " mass in 1e12 = ", Mass/1e12)

print(avg_mass_ar[M_chosen[1]]/1)
makefig.make_DM_vs_Rad_profiles_plots(series, False, dx, 999, resolution, grid_size, M_chosen, vir_rad_ar, provider, avg_mass_ar)
#makefig.make_DM_vs_Rad_profiles_plots([series[0], series[2]], True, 17, 999, resolution, grid_size, M_chosen, provider)

#make_DM_vs_Rad_profiles_plots(series, False, 15, 899)
#make_DM_vs_Rad_profiles_plots([series[0], series[2]], True, 15, 899)

#make_DM_vs_Rad_profiles_plots(series, False, 15, 799)
#make_DM_vs_Rad_profiles_plots([series[0], series[2]], True, 15, 799)
