# src/main.py

from lattice import Lattice
from monomer import Monomer
from plotter import plot_simulation

def main():
    # Step 1: Initialize the lattice and monomers
    width = 5
    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
    monomer_1 = Monomer(monomer_type='A', diffusion_rate=1.0, diffusion_energy=0.01, rotation_rate=0, rotation_energy=0, coupling_rate=1, coupling_energy=0)
    monomer_2 = Monomer(monomer_type='B', diffusion_rate=1.0, diffusion_energy=0.03, rotation_rate=0, rotation_energy=0, coupling_rate=1, coupling_energy=0)

    monomers = [monomer_1, monomer_2]

    # Step 2: Place the monomers at initial positions
    lattice.place_monomer(monomer_1, 2, 2)
    lattice.place_monomer(monomer_2, 4, 2)

    # Step 3: Call the plot_simulation function to visualize the diffusion
    plot_simulation(lattice, monomers)

if __name__ == "__main__":
    main()
