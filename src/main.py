# src/main.py

from lattice import Lattice
from monomer import Monomer
from plotter import plot_simulation, plot_final_state

def main():
    # Initialize lattice and monomers
    width = 20 # only even numbers
    monomer_params = ['A', 1.0, 0.01, 0, 0, 1, 0]

    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
    monomers = [Monomer(*monomer_params) for _ in range(50)]

    # place the monomers at initial positions
    lattice.randomly_place_monomers(monomers)

    # call the plot_simulation function to visualize the diffusion
    plot_simulation(lattice, monomers, max_steps = 100)


if __name__ == "__main__":
    main()
