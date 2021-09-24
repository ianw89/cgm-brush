from cgmbrush import *
import numpy as np


min_mass = 10**10
max_mass = 10**14.5
log_bins = 30
provider = BolshoiProvider()
load_from_files = True # If we already ran this run today, skip it
results_in_memory = False # Do not keep results in memeory, just want them saved to npy files
plots = False 
trace = False

configurations = [ 
    #Configuration(TophatProfile(), 1, provider=provider, resolution=16), 
    #Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=16), 
    #Configuration(SphericalTophatProfile(), 2, provider=provider, resolution=16), 
    #Configuration(NFWProfile(), 1, provider=provider, resolution=16),
    #Configuration(FireProfile(), 1, provider=provider, resolution=16),
    Configuration(PrecipitationProfile(), 1, provider=provider, resolution=1)
]

for config in configurations:
    config.run(trace=trace, plots=plots, load_from_files=load_from_files, results_in_memory=results_in_memory)
