# src/main.py

from lattice import Lattice
from monomer import Monomer
from plotter import plot_simulation, plot_final_state, plot_analysis_results
from analysis import analyze_structure
import random

def main():
    # Initialize lattice and monomers
    width = 100 # only even numbers

    monomer_params = ['A', 1.0, 0.00, 1.0, 0.01, 0.1, 0.000000] 
    # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy

    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True, wall=["horiz", 50, 0.0001, 0.0001])
    # monomers = [Monomer(*monomer_params) for _ in range(50)]

    # slow_growth_simulation(lattice, monomer_params, total_monomers=200, max_steps=1e6)
    wall_simulation(lattice, monomer_params, total_monomers=200, max_steps=1e6)

    # place the monomers at initial positions
    # lattice.randomly_place_monomers(monomers)

    # call the plot_simulation function to visualize the diffusion
    # plot_simulation(lattice, monomers, max_steps = 1000, animate=False)


# Approaching a more phyisically meaningful system - Reducing the complexity by adding one monomer at a time

def initialize_dimer(lattice, monomer_params):
    '''
    Here, a dimer is initialized as two coupled monomers in the centre of the lattice. The second monomer position is chosen at random 
    from the available next_nearest_neighbour positions.
    
    Returns the two monomer objects.
    '''
    x_center, y_center = (lattice.width//2, lattice.width//2) # find centre
    monomer_1 = Monomer(*monomer_params)
    monomer_1.set_position(x_center, y_center) # set position of monomer 1
    if lattice.rotational_symmetry == 6:
        orientation_1 = monomer_1.get_orientation()
        next_nearest_neighbours = lattice.get_next_nearest_neighbours(x_center, y_center, orientation_1) # get available next_nearest_neighbours positions
        monomer_2 = Monomer(*monomer_params)
        x_2, y_2 = random.choice(next_nearest_neighbours)
        monomer_2.set_position(x_2, y_2)
        monomer_2.set_orientation(orientation_1 + 180 if orientation_1 == 0 else 0) # make sure it has opposite orientation to monomer 1
        monomer_1.couple_with(monomer_2)
        lattice.place_monomer(monomer_1, x_center, y_center)
        lattice.place_monomer(monomer_2, x_2, y_2)
    return (monomer_1, monomer_2)

def introduce_new_monomer(lattice, new_monomer, monomers, first_time, max_steps=1e5):
    '''
    Introduces a new monomer on the lattice and executes the Monomer.action() method iteratively until 
    monomer is coupled (at which point it is defined to not move anymore) or until max_steps has been reached. 
    In the latter case, the monomer is removed from the lattice. This was done to avoid having several islands growing at the same time. 
    We might want to relax this condition at some point.
    '''
    steps = 0
    while steps < max_steps:
        new_monomer.action(lattice, first_time)
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
    '''
    What I call slow_growth here is what I had explained in our meeting, where the density of monomers is so low that 
    it is physically accurate to model only a single monomer at a time until it coupled to the growing island. Only after 
    the monomer has coupled is the next one introduced.
    '''
    monomer_1, monomer_2 = initialize_dimer(lattice, monomer_params)
    monomers = [monomer_1, monomer_2]
    first_time = True
    for i in range(2, total_monomers):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer]) # initialize monomer with random position (note that this can also be inside the island on an unoccupied site)
        introduce_new_monomer(lattice, new_monomer, monomers, first_time, max_steps)
        first_time = False

    print("Growth simulation completed.")

    neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)

    plot_analysis_results(neighbour_freq, radius, lattice, monomers) # some preliminary analysis of the resulting structure

def wall_simulation(lattice, monomer_params, total_monomers, max_steps=1e5):
    '''This simulation starts with a single monomer, allows it the couple to the wall, then spawns another monomer.
    '''
    monomers = []
    first_time = True
    for i in range(total_monomers):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer]) # initialize monomer with random position (note that this can also be inside the island on an unoccupied site)
        introduce_new_monomer(lattice, new_monomer, monomers, first_time, max_steps)
        first_time = False

    print("Growth simulation completed.")

    neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)

    plot_analysis_results(neighbour_freq, radius, lattice, monomers) # some preliminary analysis of the resulting structure

if __name__ == "__main__":
    main()
