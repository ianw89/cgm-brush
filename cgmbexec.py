from cgmbrush.cgmbrush import *
import numpy as np

provider = BolshoiProvider()
load_from_files = True # If we already ran this run today, skip it
results_in_memory = False # Do not keep results in memeory, just want them saved to npy files
trace = False
seed = '50g89Gfh03f4Gh0r38h2TfM08'
RS_values = RS_array_gen(1,provider.Lbox)
# TODO having all the data for all the redshifts together is too big at large resolutions
# TODO design an alternative

precip = MassDependentProfile(PrecipitationProfile(), NFWProfile(), 10**13.3)
configurations = [ 
    
    #Configuration(NFWProfile(), provider=provider, resolution=1),
    #Configuration(NFWProfile(), provider=provider, resolution=2),
    #Configuration(NFWProfile(), provider=provider, resolution=4),
    #Configuration(NFWProfile(), provider=provider, resolution=8),
    #Configuration(NFWProfile(), provider=provider, resolution=8, RS_array=RS_values),
    #Configuration(NFWProfile(), provider=provider, resolution=16),
    #Configuration(NFWProfile(), provider=provider, resolution=32),

    #Configuration(SphericalTophatProfile(), provider=provider, resolution=1),
    #Configuration(SphericalTophatProfile(), provider=provider, resolution=2),
    #Configuration(SphericalTophatProfile(), provider=provider, resolution=4),
    #Configuration(SphericalTophatProfile(), provider=provider, resolution=8),
    #Configuration(SphericalTophatProfile(), provider=provider, resolution=8, RS_array=RS_values),
    #Configuration(SphericalTophatProfile(), provider=provider, resolution=16),
    #Configuration(SphericalTophatProfile(), provider=provider, resolution=32),

    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=1),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=2),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=4),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=8),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=8, RS_array=RS_values),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=16),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=32),

    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=1),
    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=2),
    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=4),
    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=8),
    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=8, RS_array=RS_values),
    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=16),
    #Configuration(SphericalTophatProfile(rvir_factor=3), provider=provider, resolution=32),

    #Configuration(TophatProfile(), provider=provider, resolution=1),
    #Configuration(TophatProfile(), provider=provider, resolution=2),
    #Configuration(TophatProfile(), provider=provider, resolution=4),
    #Configuration(TophatProfile(), provider=provider, resolution=8),
    #Configuration(TophatProfile(), provider=provider, resolution=8, RS_array=RS_values),
    #Configuration(TophatProfile(), provider=provider, resolution=16),
    #Configuration(TophatProfile(), provider=provider, resolution=32),

    #Configuration(FireProfile(), provider=provider, resolution=1),
    #Configuration(FireProfile(), provider=provider, resolution=2),
    #Configuration(FireProfile(), provider=provider, resolution=4),
    #Configuration(FireProfile(), provider=provider, resolution=8),
    #Configuration(FireProfile(), provider=provider, resolution=8, RS_array=RS_values),
    #Configuration(FireProfile(), provider=provider, resolution=16),
    #Configuration(FireProfile(), provider=provider, resolution=32),

    Configuration(precip, provider=provider, resolution=1),
    Configuration(precip, provider=provider, resolution=2),
    Configuration(precip, provider=provider, resolution=4),
    Configuration(precip, provider=provider, resolution=8),
    Configuration(precip, provider=provider, resolution=8, RS_array=RS_values),
    Configuration(precip, provider=provider, resolution=16),
    Configuration(precip, provider=provider, resolution=32),

    #Configuration(SphericalTophatProfile(), provider=provider, resolution=32, den_grid_size=512),
    #Configuration(NFWProfile(), provider=provider, resolution=32, den_grid_size=512),
    #Configuration(FireProfile(), provider=provider, resolution=32, den_grid_size=512),
    Configuration(precip, provider=provider, resolution=32, den_grid_size=512),
    #Configuration(SphericalTophatProfile(rvir_factor=2), provider=provider, resolution=32, den_grid_size=512),

    Configuration(MassDependentProfile(SphericalTophatProfile(rvir_factor=2), NFWProfile(), 10**10.5), provider=provider, resolution=8, RS_array=RS_values),
    Configuration(MassDependentProfile(SphericalTophatProfile(rvir_factor=2), NFWProfile(), 10**11.5), provider=provider, resolution=8, RS_array=RS_values),
    Configuration(MassDependentProfile(SphericalTophatProfile(rvir_factor=2), NFWProfile(), 10**12.5), provider=provider, resolution=8, RS_array=RS_values),
    Configuration(MassDependentProfile(SphericalTophatProfile(rvir_factor=2), NFWProfile(), 10**13.5), provider=provider, resolution=8, RS_array=RS_values),
]

for config in configurations:
    config.datestamp = '2022-05-12'
    config.seed = seed

    config.run(trace=trace, load_from_files=load_from_files)
    config.generate_DM_vs_radius_profile(load_from_files=load_from_files)
    config.generate_profile_of_masks(load_from_files=load_from_files)
    
    if len(config.RS_array) > 1:
        config.generate_stacked_fields(load_from_files=load_from_files, results_in_memory=results_in_memory)
            
    # need to free up memory for next config
    config.clear_results()
