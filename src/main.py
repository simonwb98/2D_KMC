from lattice import Lattice
from monomer import Monomer

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

# Assuming Lattice and Monomer are already defined and imported

def update_grid(lattice, monomer, ax):
    """
    Updates the grid to reflect the current state of the monomer's position.
    """
    grid_data = [[0 for _ in range(lattice.width)] for _ in range(lattice.height)]

    # Mark the monomer's position with a 1 (for visualization)
    monomer_x, monomer_y = monomer.position
    grid_data[monomer_y][monomer_x] = 1

    # Clear the previous plot and redraw the grid
    ax.clear()
    ax.imshow(grid_data, cmap='Greys', vmin=0, vmax=1)
    ax.set_title(f"Monomer Position: {monomer.position}")
    ax.set_xticks(range(lattice.width))
    ax.set_yticks(range(lattice.height))
    ax.grid(True)

def run_diffusion(frame, lattice, monomer, ax):
    """
    Function to be called repeatedly by the animation to perform diffusion and update the plot.
    """
    # Diffuse the monomer to a new position
    monomer.diffuse(lattice)

    # Update the grid to show the new position
    update_grid(lattice, monomer, ax)

def main():
    # Step 1: Initialize the lattice and monomer
    width, height = 5, 5
    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
    monomer = Monomer(monomer_type='A', diffusion_rate=1.0, diffusion_energy=0.01, rotation_energy=0, rotation_rate=0)
    
    # Step 2: Place the monomer at an initial position
    lattice.place_monomer(monomer, 2, 2)

    # Step 3: Set up the plot
    fig, ax = plt.subplots()
    
    # Step 4: Use FuncAnimation to animate the diffusion process
    ani = animation.FuncAnimation(fig, run_diffusion, fargs=(lattice, monomer, ax),
                                  interval=500)  # Update every 500ms

    # Show the live plot
    plt.show()

# Run the main function to start the diffusion visualization
if __name__ == "__main__":
    main()
