# For Ian to use this on Hyak
import matplotlib.pyplot as plt
from cgmbrush import *

provider = BolshoiProvider()
STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=1, folder=varFolder)
STH_config.datestamp = '2021-09-27'
STH_config.run(load_from_files=True)
STH1_256 = STH_config.results_as_tuple

STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=2, folder=varFolder)
STH_config.datestamp = '2021-09-27'
STH_config.run(load_from_files=True)
STH2_256 = STH_config.results_as_tuple

STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=4, folder=varFolder)
STH_config.datestamp = '2021-09-27'
STH_config.run(load_from_files=True)
STH4_256 = STH_config.results_as_tuple

STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=8, folder=varFolder)
STH_config.datestamp = '2021-09-27'
STH_config.run(load_from_files=True)
STH8_256 = STH_config.results_as_tuple

STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=16, folder=varFolder)
STH_config.datestamp = '2021-09-24'
STH_config.run(load_from_files=True)
STH16_256 = STH_config.results_as_tuple

STH_config = Configuration(SphericalTophatProfile(), 1, provider=provider, resolution=32, folder=varFolder)
STH_config.datestamp = '2021-09-27'
STH_config.run(load_from_files=True)
STH32_256 = STH_config.results_as_tuple

## Image plots of methodology

# All fields:
# Halos readded field 1K, 2K, 4K, 8K


# Grid and color settings
def make_plot(vmax):

    l=0
    u=15
    du=0

    vmin = 0



    # f, ([ax1, ax2,ax4,ax5]) = plt.subplots(1, 4, figsize=(20,10)) # f is the whole plot
    #f, ([ax2,ax4,ax5,ax6]) = plt.subplots(1, 4,gridspec_kw={'hspace': 0.015, 'wspace': .01}, figsize=(20,10)) # f is the whole plot
    f, ([ax4,ax5,ax6,ax7]) = plt.subplots(1, 4,gridspec_kw={'hspace': 0.015, 'wspace': .01}, figsize=(20,10)) # f is the whole plot

    # original density

    # axis1=ax1.imshow(STH1_256[4][0,0:u,0:u],vmin=vmin,vmax=vmax)
    #axis2=ax2.imshow(STH2_256[4][0,0:2*u,0:2*u],vmin=vmin,vmax=vmax)
    axis4=ax4.imshow(STH4_256[4][0,0:4*u,0:4*u],vmin=vmin,vmax=vmax)
    axis5=ax5.imshow(STH8_256[4][0,0:8*u,0:8*u],vmin=vmin,vmax=vmax)
    axis6=ax6.imshow(STH16_256[4][0,0:16*u,0:16*u],vmin=vmin,vmax=vmax)
    axis7=ax7.imshow(STH32_256[4][0,0:32*u,0:32*u],vmin=vmin,vmax=vmax)

    # axis1=ax1.imshow(STH2_256[4][0,0:2*u,0:2*u],vmin=vmin,vmax=vmax)
    # axis2 = ax2.imshow(STH4_256[4][0,0:4*u,0:4*u],vmin=vmin,vmax=vmax)
    # axis4=ax4.imshow(fire2_256[4][0,0:2*u,0:2*u],vmin=vmin,vmax=vmax)
    # axis5=ax5.imshow(fire4_256[4][0,0:4*u,0:4*u],vmin=vmin,vmax=vmax)

    # Set ticks off
    #ax1.axis('off')
    #ax2.axis('off')
    ax4.axis('off')
    ax5.axis('off')
    ax6.axis('off')
    ax7.axis('off')


    # ax1.set_title('Resolution: 1024',size='20')
    #ax2.set_title('Resolution: 2048',size='20')
    ax4.set_title('Resolution: 4096',size='20')
    ax5.set_title('Resolution: 8192',size='20')
    ax6.set_title('Resolution: 16384',size='20')
    ax7.set_title('Resolution: 32768',size='20')

    # cbar=f.colorbar(axis1, ax=[ax1, ax2,ax4,ax5], orientation='vertical', fraction=.01,ticks=[0,vmax])
    #cbar=f.colorbar(axis2, ax=[ax2,ax4,ax5,ax6], orientation='vertical', fraction=.01,ticks=[0,vmax])
    cbar=f.colorbar(axis4, ax=[ax4,ax5,ax6,ax7], orientation='vertical', fraction=.01,ticks=[vmin,vmax])
    cbar.set_label('DM [pc cm$^{-3}$]', rotation=270,fontsize='20')

    # cbar=f.colorbar(axis1, ax=[ax1, ax2,ax4,ax5], orientation='vertical', fraction=.1,ticks=[0,vmax/3,2*vmax/3,vmax])
    # cbar.set_label('DM', rotation=0,fontsize='20')

    #f.colorbar(axis2)

    f.savefig('implot_resolution_comparison_%s.pdf' % vmax,bbox_inches='tight')


make_plot(500)
make_plot(600)
make_plot(750)
make_plot(1000)