# src/main.py

from lattice import Lattice
from monomer import Monomer
from plotter import plot_simulation, plot_final_state
import random

def main():
    # Initialize lattice and monomers
    width = 100 # only even numbers

    monomer_params = ['A', 1.0, 0.05, 1.0, 0.01, 1, 0] 
    # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy

    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
    # monomers = [Monomer(*monomer_params) for _ in range(50)]

    slow_growth_simulation(lattice, monomer_params, total_monomers=100)

    # place the monomers at initial positions
    # lattice.randomly_place_monomers(monomers)

    # call the plot_simulation function to visualize the diffusion
    # plot_simulation(lattice, monomers, max_steps = 1000, animate=False)


# Approaching a more phyisically meaningful system - Reducing the complexity by adding one monomer at a time

def initialize_dimer(lattice, monomer_params):
    x_center, y_center = (lattice.width//2, lattice.width//2)
    monomer_1 = Monomer(*monomer_params)
    monomer_1.set_position(x_center, y_center)
    if lattice.rotational_symmetry == 6:
        orientation_1 = monomer_1.get_orientation()
        next_nearest_neighbours = lattice.get_next_nearest_neighbours(x_center, y_center, orientation_1)
        monomer_2 = Monomer(*monomer_params)
        x_2, y_2 = random.choice(next_nearest_neighbours)
        monomer_2.set_position(x_2, y_2)
        monomer_2.set_orientation(orientation_1 + 180 if orientation_1 == 0 else 0)
        monomer_1.couple_with(monomer_2)
        lattice.place_monomer(monomer_1, x_center, y_center)
        lattice.place_monomer(monomer_2, x_2, y_2)
    return (monomer_1, monomer_2)

def introduce_new_monomer(lattice, new_monomer, monomers, max_steps=1e5):
    steps = 0
    while steps < max_steps:
        new_monomer.action(lattice)
        if new_monomer.coupled:
            monomers.append(new_monomer)
            print(f"Monomer succesfully coupled after {steps} steps")
            break
        steps += 1
    else:
        x, y = new_monomer.get_position()
        lattice.remove_monomer(x, y)
        print(f"Monomer failed to couple after {max_steps} steps. Initializing new monomer...")

def slow_growth_simulation(lattice, monomer_params, total_monomers, max_steps=1e5):
    monomer_1, monomer_2 = initialize_dimer(lattice, monomer_params)
    monomers = [monomer_1, monomer_2]

    for i in range(2, total_monomers):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer])
        introduce_new_monomer(lattice, new_monomer, monomers, max_steps)
        

    print("Growth simulation completed.")

    plot_final_state(lattice, monomers)

if __name__ == "__main__":
    main()
