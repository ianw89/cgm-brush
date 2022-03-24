from cgmbrush.cgmbrush import *
import numpy as np

provider = BolshoiProvider()
load_from_files = True # If we already ran this run today, skip it
results_in_memory = False # Do not keep results in memeory, just want them saved to npy files
plots = False 
trace = False
seed = '50g89Gfh03f4Gh0r38h2TfM08'
RS_values = RS_array_gen(1,provider.Lbox)
# TODO having all the data for all the redshifts together is too big at large resolutions
# TODO design an alternative

configurations = [ 
    
    Configuration(NFWProfile(), 1, provider=provider, resolution=1),
    Configuration(NFWProfile(), 1, provider=provider, resolution=2),
    Configuration(NFWProfile(), 1, provider=provider, resolution=4),
    Configuration(NFWProfile(), 1, provider=provider, resolution=8),
    Configuration(NFWProfile(), 1, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(NFWProfile(), 1, provider=provider, resolution=16),
    Configuration(NFWProfile(), 1, provider=provider, resolution=32),

    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=1),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=2),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=4),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=8),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=16),
    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=32),

    Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=32, den_grid_size=512),

    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=1),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=2),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=4),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=8),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=16),
    Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=32),

    Configuration(TophatProfile(), 1, provider=provider, resolution=1),
    Configuration(TophatProfile(), 1, provider=provider, resolution=2),
    Configuration(TophatProfile(), 1, provider=provider, resolution=4),
    Configuration(TophatProfile(), 1, provider=provider, resolution=8),
    Configuration(TophatProfile(), 1, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(TophatProfile(), 1, provider=provider, resolution=16),
    Configuration(TophatProfile(), 1, provider=provider, resolution=32),

    Configuration(FireProfile(), 1, provider=provider, resolution=1),
    Configuration(FireProfile(), 1, provider=provider, resolution=2),
    Configuration(FireProfile(), 1, provider=provider, resolution=4),
    Configuration(FireProfile(), 1, provider=provider, resolution=8),
    Configuration(FireProfile(), 1, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(FireProfile(), 1, provider=provider, resolution=16),
    Configuration(FireProfile(), 1, provider=provider, resolution=32),

    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=1),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=2),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=4),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=8),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=16),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=32),

]

for config in configurations:
    config.datestamp = '2021-11-01'
    config.seed = seed

    if len(config.RS_array) > 1:

        config.run(trace=trace, plots=plots, load_from_files=load_from_files, results_in_memory=results_in_memory)
        config.generate_DM_vs_radius_profile(load_from_files=load_from_files)
        config.generate_profile_of_masks(load_from_files=load_from_files)
        if len(config.RS_array) > 1:
            config.generate_stacked_fields(load_from_files=False, results_in_memory=results_in_memory)
                
        config.clear_results()
