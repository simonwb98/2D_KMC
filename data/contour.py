

import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

from scipy.interpolate import griddata

plt.style.use('lex_plot.mplstyle')

def plot_contour_for_experimental_file(results_csv, experimental_file):

    """Generate contour plots for each rotation energy for a given experimental file."""

    # Load the results CSV

    results_data = pd.read_csv(results_csv)


    # Filter the data for the specific experimental file

    file_data = results_data[results_data["Experimental File"] == experimental_file]


    # Get unique Rotation Energy values

    rotation_energies = file_data["Dehalogenation Energy"].unique()
    print(rotation_energies)
    i = 0
    fig = plt.figure()
    
#     ax1 = fig.add_axes((.1-0.05,.1, .8/3, .8))
#     #ax1.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
#     
#     ax1.set_ylabel('Coupling Energy [eV]')
#     #ax1.yaxis.set_label_coords(-.1, .0)
#     #ax1.text(.04, 1.9, r"$W_d$",fontsize='large')
#     ax1.set_xticks(np.linspace(0.55, 0.95, 5))
#     
#     ax2 = fig.add_axes((.1 +.8/3-0.05, .1, .8/3, .8))
#     ax2.set_yticks([])
#     ax2.set_xticks(np.linspace(0.55, 0.95, 5))
#     ax2.set_xlabel('Diffusion Energy [eV]')
#     
#     ax3 = fig.add_axes((.1 +2*.8/3-0.05, .1, .8/3+0.033, .8))
#     ax3.set_yticks([])
#     ax3.set_xticks(np.linspace(0.55, 0.95, 5))
    
    
    ax1 = fig.add_axes((.1-0.05,.5, .8/3, .4))
    ax1.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
    ax1.set_xticks([])
    ax1.set_ylabel('Coupling Energy [eV]')
    ax1.yaxis.set_label_coords(-.1, .0)
    #ax1.text(.04, 1.9, r"$W_d$",fontsize='large')
    
    ax2 = fig.add_axes((.1 +.8/3-0.05, .5, .8/3, .4))
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax3 = fig.add_axes((.1 +2*.8/3-0.05, .5, .8/3+0.033, .4))
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax4 = fig.add_axes((.1-0.05,.1, .8/3, .4))
    ax4.set_yticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
    ax4.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
    #ax4.text(.04, -0.5, r"$W_d$",fontsize='large')
    
    ax5 = fig.add_axes((.1 +.8/3-0.05,.1, .8/3, .4))
    ax5.set_yticks([])
    ax5.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
    ax5.set_xlabel('Diffusion Energy [eV]')
    
    ax6 = fig.add_axes((.1 +2*.8/3 -0.05,.1, .8/3+0.033, .4))
    ax6.set_yticks([])
    ax6.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]

    #axes = [ax1, ax2, ax3]
    levels = np.linspace(0.0, 0.4, 21)
    
    for rotation_energy in rotation_energies:
        
        # Filter data by Rotation Energy

        rotation_data = file_data[file_data["Dehalogenation Energy"] == rotation_energy]
        

        # Extract Diffusion Energy, Coupling Energy, and Wasserstein Distance

        x = rotation_data["Diffusion Energy"].values

        y = rotation_data["Coupling Energy"].values

        z = rotation_data["Wasserstein Distance"].values

        
        # Create a grid for contour plotting

        x_grid = np.linspace(min(x), max(x), 100)

        y_grid = np.linspace(min(y), max(y), 100)

        x_grid, y_grid = np.meshgrid(x_grid, y_grid)
        

        # Interpolate z values onto the grid

        z_grid = griddata((x, y), z, (x_grid, y_grid), method="cubic")
        print(min(z_grid.reshape(100**2)))
        print(np.where(z_grid == min(z_grid.reshape(100**2))))
        
        # Generate the contour plot

        #plt.figure(figsize=(10, 8))
        
        contour = axes[i].contourf(x_grid, y_grid, z_grid, cmap="Oranges", levels=levels)
        axes[i].text(.1, 1.55, r'$E_{De}=$'+ str(rotation_energy) + 'eV')
#         if i in [0, 2]:
#             cont = fig.colorbar(contour, location='top',ax=axes[i],shrink=0.7)
#             cont.set_ticks(ticks = np.linspace(min(z), max(z), 6), labels=np.around(np.linspace(min(z), max(z), 6),2))
#         elif i == 1:
#             cont = fig.colorbar(contour, location='top',ax=axes[i],shrink=0.7)
#             cont.set_ticks(ticks = np.linspace(min(z), max(z), 6), labels=np.around(np.linspace(min(z), max(z), 6),2))
#             
#         else:
#             cont = fig.colorbar(contour, location='bottom',ax=axes[i],shrink=0.7)
#             cont.set_ticks(ticks = np.linspace(min(z), max(z), 6), labels=np.around(np.linspace(min(z), max(z), 6),2))
#         
        if i == 2:
            cont = fig.colorbar(contour, label='Wasserstein Distance', location='right',ax=[axes[2], axes[-1]],\
                                shrink=1, pad=0.0)
            #cont.set_ticks(ticks = np.linspace(0.03, 0.3, 7),labels=np.around(np.linspace(0.03, 0.3, 7),2))
        #plt.xlabel("Diffusion Energy")

        #plt.ylabel("Coupling Energy")

        #plt.title(f"Contour Plot for {experimental_file}\nDehalogenation Energy = {rotation_energy}")

        #plt.show()
        i += 1