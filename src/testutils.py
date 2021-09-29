from cgmbrush import *
import numpy as np
import matplotlib.pyplot as plt
#from truth.truth import AssertThat # considering using PyTruth for fluent assertions

def plot_grid_comparison(new, baseline):
    fig, ax = plt.subplots(1,2,figsize=(24, 12))

    vmin = 0
    vmax = max(np.max(new), np.max(baseline))
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    pos = ax[0].imshow(new, norm=norm)
    ax[0].title.set_text('New Result')
    pos = ax[1].imshow(baseline, norm=norm)
    ax[1].title.set_text('Baseline')

    fig.colorbar(pos, ax=ax, shrink=0.85)

    return fig

def compare_mask_lists(results, baseline, resolution):
    for i in range(len(results)):
        if not np.allclose(results[i], baseline[i]):
            print("Mask %s changed." % i)
            plot_grid_comparison(results[i], baseline[i])
            print("Values of center pixel:")
            print(" New:\t\t %s" % results[i][10*resolution][10*resolution])
            print(" Baseline:\t %s" % baseline[i][10*resolution][10*resolution])
            