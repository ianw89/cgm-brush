from cgmbrush.cgmbrush import *
from cgmbrush.plots.plots import *
import matplotlib.pyplot as plt
import numpy as np

provider = BolshoiProvider()

resolutions = [4,8,16,32]
profiles = [SphericalTophatProfile(), NFWProfile(), FireProfile(), PrecipitationProfile()]

configs = []

for p in profiles:
    for r in resolutions:
        configs.append(Configuration(p, 1, provider=provider, resolution=r))

summary = ''

for c in configs:
    c.datestamp = '2022-04-04' 
    c.run(load_from_files=True)
    std = DM_statistics(c.get_final_field()[0])
    summary = summary + "{} ({}):\t{}\n".format(c.addition_profile.pretty_name, c.resolution, std)
    c.clear_results()

print(summary)