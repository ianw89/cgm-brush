###################################################################################################
#
# plots.py 	                (c) Ian Williams, Adnan Khan, Matt McQuinn
#     				    	ianw89@live.com
#
###################################################################################################

from cgmbrush.cgmbrush import *
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Rectangle

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 16
XBIG_SIZE = 20

# Plotting Defaults
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=XBIG_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def view_mask(profile: CGMProfile, mass):

    provider = BolshoiProvider()

    # pick some parameters
    redshift = 0 
    resolution = 16
    cellsize = provider.Lbox / (provider.halofieldresolution * resolution)
    #print("cellsize (kpc): %s" % (cellsize * 1000))
    fine_mask_len = 20 * resolution
    zs = fine_mask_len - fine_mask_len // 10
    ze = fine_mask_len + fine_mask_len // 10

    comoving_rvir = halo.comoving_rvir(cosmo, mass, redshift)
    #print("co-moving virial radius (kpc): %s" % (comoving_rvir * 1000))

    mask, area = normalize_and_convert_mask(profile.get_mask(mass, comoving_rvir, redshift, resolution, cellsize, fine_mask_len), mass, cellsize)
    fig, ax = plt.subplots(1,2,figsize=(24, 12))
    pos = ax[0].imshow(mask)
    pos = ax[1].imshow(mask[zs:ze,zs:ze])
    ax[0].title.set_text('{} Profile, M=${:.1e}$'.format(profile.name, mass))
    ax[1].title.set_text('{} Profile (zoomed), M=${:.1e}$'.format(profile.name, mass))
    cbar = fig.colorbar(pos, ax=ax)
    cbar.set_label('DM [pc cm$^{-3}$]', rotation=270, size=XBIG_SIZE)


def density_plot(axis, field, vmin: int, vmax: int, name: str):
    a = axis.imshow(field, vmin = vmin, vmax = vmax)
    axis.set_title(name, size='20')
    return a

def fields_comparison_plot(fields, vmin, vmax):
    """Creates a plot 1 x n long for n fields with a single shared color bar."""
    n = len(fields)

    fig,axs = plt.subplots(1, n, gridspec_kw={'hspace': 0.015, 'wspace': .01}, figsize=(5*n,5))
    
    for i in range(0, n):
        field = fields[i][1]
        name = fields[i][0]
        a1 = density_plot(axs[i], field, vmin, vmax, name)
        axs[i].axis('off')

    cbar = fig.colorbar(a1, ax=axs, orientation='vertical', fraction=.05,ticks=[vmin,vmax])
    cbar.set_label('DM [pc cm$^{-3}$]', rotation=270, fontsize='20')


    return fig,axs


def radial_distances_kpc(cellsize, num_values):
    MpctoKpc =1000
    return MpctoKpc * (0.5 + np.arange(0, num_values)) * cellsize    


def make_DM_vs_Rad_profiles_plots(series, plot_masks: bool, x_start, x_end, resolution, grid_size, M_chosen,  vir_rad_ar, provider, avg_mass_ar, ylims, name, y_label='DM - <DM> [pc cm$^{-3}$]', customFunc=None):
    ''' This plots the radial profile.  
    error = True makes same assumptions as Khan++2022 to compute the error  (data from a precomputation of the error of PDF is included)
    xstart and xend
    M_chosen are mass bins
    '''

    mean_DM = np.mean(provider.get_density_field(0, 256, 2))

    plots_to_make = len(M_chosen)

    XXBIG_SIZE = 30
    MASSIVE_SIZE = 50

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=XXBIG_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=XXBIG_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=XBIG_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=XXBIG_SIZE)  # fontsize of the figure title

    MpctoKpc = 1000
    hspace = 0.015

    DM_Rad_fig, DM_Rad_axs = plt.subplots(plots_to_make, 1,
                            gridspec_kw={'hspace': hspace, 'wspace': .2},figsize=(20,plots_to_make*10))

    if plots_to_make == 1:
        DM_Rad_axs = [DM_Rad_axs]

    dx = (provider.Lbox/grid_size)  #cell size
    x_axis = radial_distances_kpc(dx, series[0][0].shape[1])

    for i in range(plots_to_make):
        plot_DM_vs_Rad(x_axis, mean_DM, DM_Rad_axs[i], M_chosen[i], series, plot_masks)

    # Hook for callers to add custom plotting code
    if customFunc is not None:
        customFunc(series, M_chosen, DM_Rad_axs)
    
    # ticks

    # ax.xaxis.grid(True, which='minor')
    # DM_Rad_axs[0].xaxis.grid(axis='x', which='minor', bottom=True)

    # DM_Rad_axs[0].tick_params(axis='x', which='minor', bottom=True)
    # DM_Rad_axs[1].tick_params(axis='x', which='minor', bottom=True)

    # DM_Rad_axs[0].tick_params(axis='x', which='minor',fontsize=10)
    # DM_Rad_axs[2].xaxis.set_tick_params(width=5)

    for i in range(plots_to_make):
        DM_Rad_axs[i].tick_params('both', length=10, width=4, which='major')
        DM_Rad_axs[i].tick_params('both', length=10, width=4, which='minor')
        # x and y axis limits
        DM_Rad_axs[i].set_xlim([x_start, x_end])
        DM_Rad_axs[i].set_ylim([0, ylims[i]])

    # legend
    DM_Rad_axs[0].legend(loc='right',prop={'size':28}, frameon=False)
    DM_Rad_axs[math.floor((plots_to_make-1)/2)].set_ylabel(y_label,fontsize=MASSIVE_SIZE)
    DM_Rad_axs[plots_to_make-1].set_xlabel('Impact Parameter [kpc]',fontsize=MASSIVE_SIZE)

    # Adding virial radii vertical line
    if (vir_rad_ar is not None):
        for i in range(plots_to_make):
            DM_Rad_axs[i].axvline(1*MpctoKpc*(vir_rad_ar[M_chosen[i]]), color='k', linestyle='--', linewidth=1)
            #DM_Rad_axs[i].axvline(2*MpctoKpc*(vir_rad_ar[M_chosen[i]]), color='k', linestyle='--', linewidth=1)
            #DM_Rad_axs[i].axvline(3*MpctoKpc*(vir_rad_ar[M_chosen[i]]), color='k', linestyle='--', linewidth=1)
            
    # Rectangular patch for region too far inside resolution limit
    #DM_Rad_axs[0].add_patch(Rectangle((0,0), 45, 130,facecolor='yellow'))
    #DM_Rad_axs[1].add_patch(Rectangle((0,0), 45, 399,facecolor='yellow'))
    #DM_Rad_axs[2].add_patch(Rectangle((0,0), 45, 1500,facecolor='yellow'))

    # mass labels
    for i in range(plots_to_make):
        DM_Rad_axs[i].text(x_end*0.4, ylims[i]*0.88, r'$10^{%.1f} M_\odot$ '%np.log10(avg_mass_ar[M_chosen[i]]),fontsize=34) #'Mass = %.1E $M_\odot$' % Decimal(df[2][M_chosen[1]]),fontsize=30)

    # DM_Rad_axs[0].rc('xtick', labelsize=35)    # fontsize of the tick labels
    # # plt.rc('ytick', labelsize=35)    # fontsize of the tick labels

    name = 'DMvsRad_{}_{}_{}.pdf'.format(name, resolution, x_end)

    saveFig(name, DM_Rad_fig, bbox_inches='tight')

def plot_DM_vs_Rad(x_axis, mean_DM, axis, massbin, series, plot_masks):
    
    # axis.set_title('Mass = %.1E' % Decimal(df[2][massbin]),fontsize=14)
    
    for i in range(len(series)):
        data = series[i]
        #print("terms = ", mean_DM, data[0][massbin,:]-mean_DM, x_axis,data[1][massbin,:])
        axis.semilogx(x_axis,data[0][massbin,:]-mean_DM, data[4], label=data[2],lw=5, color=data[3])

        if plot_masks:
            axis.semilogx(x_axis,data[1][massbin,:], '--', lw=2,color=data[3])


def plot_error_bars(r_star_ar, axis, *lines):
    for line in lines:
        axis.semilogx(r_star_ar*1000,line[0],ls='--',drawstyle='steps',color=line[1],lw=line[2])

    

def error_bar_DMvsRad(sd,N_frb,radius_array,Rvir,avg_frbs):
    """Function that computes error bar for plot
    avg_frbs is how many halos of this size an FRB goes through; see McQuinn 2014."""
    error_bar = np.zeros([len(radius_array)])
    for i in range(1,len(radius_array)):
        error_bar[i-1] = sd/(np.sqrt(avg_frbs*N_frb*((radius_array[i]**2-radius_array[i-1]**2)/Rvir**2)))
    
    return error_bar

