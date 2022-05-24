from cgmbrush.cgmbrush import *
from cgmbrush.plots.plots import *
import matplotlib.pyplot as plt
import numpy as np

provider = BolshoiProvider()
date = '2022-04-04'
date2 = '2022-05-12' 
resolutions = [4,8,16,32]
precipProfile = MassDependentProfile(PrecipitationProfile(), NFWProfile(), 10**13.3)
profiles = [(SphericalTophatProfile(), date), (SphericalTophatProfile(rvir_factor=2), date), (NFWProfile(), date), (FireProfile(), date2), (precipProfile, date2)]

configs = []

for p in profiles:
    for r in resolutions:
        configs.append(Configuration(p[0], provider=provider, resolution=r, datestamp=p[1]))

summary = ''

configs[0]

for c in configs:
    c.datestamp = '2022-05-12' 
    c.run(load_from_files=True)
    std = DM_statistics(c.get_final_field()[0])
    summary = summary + "{} ({}):\t{}\n".format(c.addition_profile.pretty_name, c.resolution, std)
    c.clear_results()

print(summary)