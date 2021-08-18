from cgmbrush import *
import numpy as np


min_mass = 10**10
max_mass = 10**14.5
log_bins = 30
config = Configuration('tophat_spherical', 1, resolution=4, file_prefix='STH', load_from_files=False)
config.provider = BolshoiProvider()

config.results = hist_profile(config.provider, config.den_grid_size, config.RS_array, min_mass, 
                                    max_mass, log_bins, config.subtraction_halo_profile, 
                                    config.addition_halo_profile, config.scaling_radius, config.resolution)

original = config.provider.get_density_field(0, config.den_grid_size)
background_dm = config.results[5][0]
cgm_only = config.results[4][0]
density_final = config.results[1][0]

fig, axes = plt.subplots(1,4,figsize=(18, 4))
pos = axes[0].imshow(original) 
fig.colorbar(pos, ax=axes[0])
axes[0].title.set_text('Original Density Field')
pos = axes[1].imshow(background_dm) 
fig.colorbar(pos, ax=axes[1])
axes[1].title.set_text('Density minus halos')
pos = axes[2].imshow(cgm_only) 
fig.colorbar(pos, ax=axes[2])
axes[2].title.set_text('CGM Profile to add')
pos = axes[3].imshow(density_final) 
fig.colorbar(pos, ax=axes[3])
axes[3].title.set_text('Final Product')

filename = config.file_prefix + str(config.resolution) + '_' + str(config.den_grid_size)
saveFig(filename, fig)
saveArray(filename, *config.results)