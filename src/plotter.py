# src/plotter.py

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import RegularPolygon, Circle
import math

def update_hexagonal_grid(lattice, monomers, ax):
    """
    Updates the hexagonal grid to reflect the current state of the monomer's position.
    Each hexagonal grid point is represented by a circle.
    """
    ax.clear()

    # Loop through all valid lattice coordinates and plot circles at each point
    for (x, y) in lattice.lattice_coord:
        # Adjust x-position to create the staggered hexagonal effect
        x_offset = x + 0.5 if y % 2 != 0 else x

        ax.add_patch(plt.Circle((x_offset, y), 0.3, facecolor='lightgray', edgecolor='black', lw=1))

    # Plot the monomers' positions
    for monomer in monomers:
        x_mon, y_mon = monomer.position
        x_mon = x_mon + 0.5 if y_mon % 2 != 0 else x_mon

        orientation = monomer.get_orientation()

        color = "blue" if monomer.coupled else "red"
        triangle = RegularPolygon(
            (x_mon, y_mon),
            numVertices=3,
            radius=0.6,
            orientation=math.radians(orientation),
            facecolor=color,
            edgecolor="black",
            lw=2
        )
        
        ax.add_patch(triangle)

    ax.set_xlim(-0.5, lattice.width)
    ax.set_ylim(-0.5, lattice.height)
    ax.set_aspect('equal')  # Keep the aspect ratio so circles don't look squashed
    ax.grid(False)  # Disable the grid lines
    ax.set_xticks([])
    ax.set_yticks([])

def run_actions(lattice, monomers):
    """Perform diffusion and update the hexagonal plot."""
    for monomer in monomers:
        monomer.action(lattice)

def update_plot(frame, lattice, monomers, ax):
    run_actions(lattice, monomers)
    update_hexagonal_grid(lattice, monomers, ax)

def all_monomers_coupled(monomers):
    # stopping condition - when all monomers have coupled together
    return all(monomer.coupled for monomer in monomers)

def plot_simulation(lattice, monomers, max_steps=1000, animate=False):
    """
    Sets up the plot and runs the animation based on the current state of the simulation.
    """
    steps = 0
    while steps < max_steps:
        run_actions(lattice, monomers)

        if all_monomers_coupled(monomers):
            print(f"Termination condition reached after {steps + 1} steps.")
            break

        steps += 1

    if animate:
        fig, ax = plt.subplots()

        ani = animation.FuncAnimation(fig, update_plot, fargs=(lattice, monomers, ax), interval=10, frames = 20) #  blit=True

        plt.show()
    else:
        plot_final_state(lattice, monomers)

def plot_final_state(lattice, monomers):
    fig, ax = plt.subplots()
    update_hexagonal_grid(lattice, monomers, ax)
    plt.show()

def plot_analysis_results(neighbor_frequencies, effective_radius, lattice, monomers, defects=False, herringbone=False):
    """Plot the polymer structure with the radius of gyration and the histogram of coupled neighbors."""
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 7)) 

    ax[0].clear()
    
    # loop through all valid lattice coordinates and plot circles at each point
    for (x, y) in lattice.lattice_coord:
        x_offset = x + 0.5 if y % 2 != 0 else x  # adjust for hexagonal staggered rows
        ax[0].add_patch(plt.Circle((x_offset, y), 0.3, facecolor='lightgray', edgecolor='black', lw=1))

    # plot the monomers' positions
    for monomer in monomers:
        x_mon, y_mon = monomer.position
        if y_mon % 2 != 0:
            x_mon += 0.5  # offset for odd rows

        orientation = monomer.get_orientation()
        color = "blue" if monomer.coupled else "red"
        triangle = RegularPolygon(
            (x_mon, y_mon),
            numVertices=3,
            radius=0.7,
            orientation=math.radians(orientation),
            facecolor=color,
            edgecolor="black",
            lw=2
        )
        ax[0].add_patch(triangle)
    
    # plot the defects' positions
    if defects:
        for defect in defects:
            x, y = defect.position
            if y % 2 != 0:
                x += 0.5  # offset for odd rows

            color = 'green' if defect.nucleating else 'red'
            
            circ = RegularPolygon(
                (x, y),
                numVertices=8,
                radius=0.3,
                facecolor=color,
                edgecolor="black",
                lw=2
            )
            ax[0].add_patch(circ)
        
    if herringbone:
        for defect in herringbone:
            x, y = defect.position
            if y % 2 != 0:
                x += 0.5  # offset for odd rows

            color = 'green' if defect.nucleating else 'red'
            
            circ = RegularPolygon(
                (x, y),
                numVertices=8,
                radius=0.3,
                facecolor=color,
                edgecolor="black",
                lw=2
            )
            ax[0].add_patch(circ)
        
    # calculate the center of mass for the polymer
    x_positions = [monomer.position[0] for monomer in monomers]
    y_positions = [monomer.position[1] for monomer in monomers]
    
    x_center_of_mass = sum(x_positions) / len(x_positions)
    y_center_of_mass = sum(y_positions) / len(y_positions)

    # add the circle representing the radius of gyration
    circle = Circle((x_center_of_mass, y_center_of_mass), effective_radius, fill=False, color='green', lw=2, linestyle='--')
    ax[0].add_patch(circle)
    
    ax[0].set_title('Polymer Structure with Radius of Gyration')
    ax[0].set_xlim(-0.5, lattice.width)
    ax[0].set_ylim(-0.5, lattice.height)
    ax[0].set_aspect('equal')
    ax[0].grid(False)
    ax[0].set_xticks([]) # get rid of tick marks
    ax[0].set_yticks([])

    # plot the histogram of neighbor frequencies on the second axis
    labels = sorted(neighbor_frequencies.keys())
    counts = [neighbor_frequencies[label] for label in labels]
    
    ax[1].bar(labels, counts, color='blue', edgecolor='black')
    ax[1].set_title('Histogram of Coupled Neighbors')
    ax[1].set_xlabel('Number of Coupled Neighbors')
    ax[1].set_ylabel('Monomer Count')

    plt.tight_layout()
    plt.show()
