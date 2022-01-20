###################################################################################################
#
# testutils.py 	        (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################

from cgmbrush.cgmbrush import *
import numpy as np
import matplotlib.pyplot as plt
#from truth.truth import AssertThat # considering using PyTruth for fluent assertions

def plot_grid_comparison(new, baseline):
    fig, ax = plt.subplots(1,2,figsize=(24, 12))

    if np.shape(new) is not np.shape(baseline):
        print("Differing shapes. Baseline: {0}.  New: {1}.".format(np.shape(baseline), np.shape(new)))

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
    fig = None
    for i in range(len(results)):
        if not np.allclose(results[i], baseline[i]):
            print("Mask %s changed." % i)
            print("Values of center pixel:")
            print(" New:\t\t %s" % results[i][10*resolution][10*resolution])
            print(" Baseline:\t %s" % baseline[i][10*resolution][10*resolution])
            fig = plot_grid_comparison(results[i], baseline[i])
    
    return fig
    

def check_validity(config: Configuration):
    #assert not np.isnan(results['massbin_histograms']), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_subtraction_coarse_field())), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_removed_field())), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_addition_masks())), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_addition_field())), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_final_field())), 'No NaN should appear in results.'
    #assert not np.any(np.isnan(results['stacked_density_field'])), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_virial_radii())), 'No NaN should appear in results.'
    assert not np.any(np.isnan(config.get_halo_masses())), 'No NaN should appear in results.'

def check_stacked_validity(config: Configuration):
    assert not np.any(np.isnan(config.get_stacked_orig_field())), 'No NaN should appear in results.'
    assert not np.allclose(config.get_stacked_orig_field(), 0), 'Fields shouldn\'t just be 0\'s.'
    assert not np.any(np.isnan(config.get_stacked_removed_field())), 'No NaN should appear in results.'
    assert not np.allclose(config.get_stacked_removed_field(), 0), 'Fields shouldn\'t just be 0\'s.'
    assert not np.any(np.isnan(config.get_stacked_addition_field())), 'No NaN should appear in results.'
    assert not np.allclose(config.get_stacked_addition_field(), 0), 'Fields shouldn\'t just be 0\'s.'
    assert not np.any(np.isnan(config.get_stacked_final_field())), 'No NaN should appear in results.'
    assert not np.allclose(config.get_stacked_final_field(), 0), 'Fields shouldn\'t just be 0\'s.'



def force_load_npz(filename, folder=VAR_DIR):
    """Loads all numpy arrays in a npz into a dictionary."""
    file_path = os.path.join(folder, filename + ".npz")
    results = {}
    npz = np.load(file_path, allow_pickle=True)
    for file in npz.files:
        results[file] = npz[file] 
    npz.close()
    return results