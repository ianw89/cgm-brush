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

def algorithm_stages_plot(raw_field, removed_field, addition_field, final_field, orig_grid_size, zoom_rect: bool):
    vmin = 0
    vmax = 300
    fig, axs = fields_comparison_plot( [('Raw density field', raw_field), ('Halos removed', removed_field), ('Halo-associated gas', addition_field), ('Final halo density field', final_field)], vmin, vmax)
    
    if zoom_rect:
        rect = Rectangle((0,0),320,320,linewidth=1,edgecolor='r',facecolor="none")
        axs[3].add_patch(rect)
    
    filename = 'implot_add_remove_halos_{}.pdf'.format(orig_grid_size)
    saveFig(filename, fig, bbox_inches='tight')



def compare_4_profile_fields(STH2, NFW, FIRE, PRE, resolution, z):
    """Creates a plot comparing the fields from 4 CGM profiles (see variable names)."""

    vmin = 0
    vmax = max(np.max(STH2), np.max(NFW), np.max(FIRE), np.max(PRE))
    vhigh = vmax*0.25
    vmax = int(round(vhigh, -int(math.floor(math.log10(vhigh)))))  


    fig, axs = fields_comparison_plot([('NFW', NFW), ('FIRE', FIRE), ('Precipitation', PRE), ('Tophat $2R_{vir}$', STH2)], vmin, vmax)

    filename = 'stacked_field_images_{}_z{:.1f}.pdf'.format(resolution, z)
    saveFig(filename, fig, bbox_inches='tight')
