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

    if np.shape(new) != np.shape(baseline):
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

def plot_grid_comparison_delta(new, baseline):

    fig, ax = plt.subplots(1,3,figsize=(36, 12))

    if np.shape(new) != np.shape(baseline):
        print("Differing shapes. Baseline: {0}.  New: {1}.".format(np.shape(baseline), np.shape(new)))

    vmin = 0
    vmax = max(np.max(new), np.max(baseline))
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    pos = ax[0].imshow(new, norm=norm)
    ax[0].title.set_text('New Result')
    fig.colorbar(pos, ax=ax[0])
    pos = ax[1].imshow(baseline, norm=norm)
    ax[1].title.set_text('Baseline')
    fig.colorbar(pos, ax=ax[1])

    pos = ax[2].imshow(new - baseline) # delta between them
    fig.colorbar(pos, ax=ax[2])

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

def add_cross(field, center_x, center_y, weight, falloff=3):
    field[center_x][center_y] += weight
    field[center_x - 1][center_y] += weight/falloff
    field[center_x + 1][center_y] += weight/falloff
    field[center_x][center_y - 1] += weight/falloff
    field[center_x][center_y + 1] += weight/falloff

def add_big_cross(field, center_x, center_y, weight):
    field[center_x][center_y] += weight
    field[center_x - 1][center_y] += weight/3
    field[center_x + 1][center_y] += weight/3
    field[center_x][center_y - 1] += weight/3
    field[center_x][center_y + 1] += weight/3
    field[center_x - 2][center_y] += weight/10
    field[center_x + 2][center_y] += weight/10
    field[center_x][center_y - 2] += weight/10
    field[center_x][center_y + 2] += weight/10
    field[center_x + 1][center_y + 1] += weight/10
    field[center_x - 1][center_y + 1] += weight/10
    field[center_x + 1][center_y - 1] += weight/10
    field[center_x - 1][center_y - 1] += weight/10

class FakeProvider(SimulationProvider): 

    Lbox = 50 / cosmo.h # 50 Mpc/h fake box
    halofieldresolution = 50 # unlike bolshoi making this match the original density field size
    
    def get_density_field(self, redshift: float, resolution: int, proj_axis: int):
        return np.zeros((50,50)) # 50 x 50 cell original grid

    def get_halos(self, redshift : int) -> pd.DataFrame:
        d = {'x': [3,15,13,30,30], 'y': [3,10,13,37,40], 'Mvir': [10**12, 10**12, 10**12, 10**12, 10**13]} # x,y coords are in Mpc/h like Bolshoi
        df = pd.DataFrame(data=d)
        return df

    def get_z_name(redshift: float) -> str:
        """Gets the shorthand name associated with a given redshift, used in filenames."""
        raise NotImplementedError