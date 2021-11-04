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

