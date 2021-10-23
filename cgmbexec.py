from cgmbrush.cgmbrush import *
import numpy as np


min_mass = 10**10
max_mass = 10**14.5
log_bins = 30
provider = BolshoiProvider()
load_from_files = True # If we already ran this run today, skip it
results_in_memory = False # Do not keep results in memeory, just want them saved to npy files
plots = False 
trace = False
RS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# TODO having all the data for all the redshifts together is too big at large resolutions
# TODO design an alternative

configurations = [ 
    Configuration(NFWProfile(), 1, provider=provider, resolution=1),
    Configuration(NFWProfile(), 1, provider=provider, resolution=2),
    Configuration(NFWProfile(), 1, provider=provider, resolution=4),
    Configuration(NFWProfile(), 1, provider=provider, resolution=8),
    Configuration(NFWProfile(), 1, provider=provider, resolution=8, RS_array=RS),
    Configuration(NFWProfile(), 1, provider=provider, resolution=16),
    Configuration(NFWProfile(), 1, provider=provider, resolution=32),

    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=1),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=2),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=4),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=8),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=8, RS_array=RS),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=16),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=32),

    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=1),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=2),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=4),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=8),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=8, RS_array=RS),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=16),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=32),

    Configuration(TophatProfile(), 1, provider=provider, resolution=1),
    Configuration(TophatProfile(), 1, provider=provider, resolution=2),
    Configuration(TophatProfile(), 1, provider=provider, resolution=4),
    Configuration(TophatProfile(), 1, provider=provider, resolution=8),
    Configuration(TophatProfile(), 1, provider=provider, resolution=8, RS_array=RS),
    Configuration(TophatProfile(), 1, provider=provider, resolution=16),
    Configuration(TophatProfile(), 1, provider=provider, resolution=32),

    Configuration(FireProfile(), 1, provider=provider, resolution=1),
    Configuration(FireProfile(), 1, provider=provider, resolution=2),
    Configuration(FireProfile(), 1, provider=provider, resolution=4),
    Configuration(FireProfile(), 1, provider=provider, resolution=8),
    Configuration(FireProfile(), 1, provider=provider, resolution=8, RS_array=RS),
    Configuration(FireProfile(), 1, provider=provider, resolution=16),
    Configuration(FireProfile(), 1, provider=provider, resolution=32),

    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=1),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=2),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=4),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=8),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=8, RS_array=RS),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=16),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=32),

]

for config in configurations:
    config.run(trace=trace, plots=plots, load_from_files=load_from_files, results_in_memory=results_in_memory)
    config.generate_DM_vs_radius_profile()
    config.generate_profile_of_masks()
    config.generate_stacked_field(results_in_memory=results_in_memory)
