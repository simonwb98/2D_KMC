# src/plotter.py

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import RegularPolygon
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
        if y % 2 != 0:
            x_offset = x + 0.5
        else:
            x_offset = x

        ax.add_patch(plt.Circle((x_offset, y), 0.3, facecolor='lightgray', edgecolor='black', lw=1))

    # Plot the monomers' positions
    for monomer in monomers:
        x_mon, y_mon = monomer.position
        if y_mon % 2 != 0:
            x_mon += 0.5  # Apply offset for odd rows

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

def all_monomers_coupled(monomers):
    # stopping condition when all monomers have coupled together
    return all(monomer.coupled for monomer in monomers)

def plot_simulation(lattice, monomers, max_steps=1000, animation=False):
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

    if animation:
        fig, ax = plt.subplots()

        # Use FuncAnimation to animate the diffusion process
        ani = animation.FuncAnimation(fig, run_actions, fargs=(lattice, monomers, ax), interval=10, frames = 20, save_count=10)

        # Show the live plot
        plt.show()
    else:
        plot_final_state(lattice, monomers)

def plot_final_state(lattice, monomers):
    fig, ax = plt.subplots()
    update_hexagonal_grid(lattice, monomers, ax)
    plt.show()