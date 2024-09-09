import matplotlib.pyplot as plt
import matplotlib.animation as animation

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
        monomer_x, monomer_y = monomer.position
        if monomer_y % 2 != 0:
            monomer_x += 0.5  # Apply offset for odd rows

        ax.add_patch(plt.Circle((monomer_x, monomer_y), 0.3, facecolor = 'red', edgecolor='black', lw=2))

    ax.set_xlim(-0.5, lattice.width)
    ax.set_ylim(-0.5, lattice.height)
    ax.set_aspect('equal')  # Keep the aspect ratio so circles don't look squashed
    ax.grid(False)  # Disable the grid lines
    ax.set_xticks([])
    ax.set_yticks([])

def run_diffusion(frame, lattice, monomers, ax):
    """Perform diffusion and update the hexagonal plot."""
    for monomer in monomers:
        monomer.diffuse(lattice)
    update_hexagonal_grid(lattice, monomers, ax)

def plot_simulation(lattice, monomers):
    """
    Sets up the plot and runs the animation based on the current state of the simulation.
    """
    fig, ax = plt.subplots()

    # Use FuncAnimation to animate the diffusion process
    ani = animation.FuncAnimation(fig, run_diffusion, fargs=(lattice, monomers, ax), interval=100, save_count=100)

    # Show the live plot
    plt.show()
