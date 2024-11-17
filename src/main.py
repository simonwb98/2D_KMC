# src/main.py

from lattice import Lattice
from monomer import Monomer
from plotter import plot_simulation, plot_final_state, plot_analysis_results
from analysis import analyze_structure
import random
import numpy as np
from matplotlib import pyplot as plt

def main():
    N = 2 # N^2 is the number of simulations

    p_c = np.linspace(0.1, 1, N)
    p_d = np.linspace(0.1, 1, N)

    A, B = np.meshgrid(p_c, p_d)
    save_meshgrid(A, B, 
                  r"C:\Users\User\Desktop\research-updates\thesis\figures\results\2D_KMC\meshgrid_p_c.npy", 
                  r"C:\Users\User\Desktop\research-updates\thesis\figures\results\2D_KMC\meshgrid_p_d.npy")
    print("Meshgrids p_c and p_d saved.")

    results = np.zeros((N, N), dtype=[
        ('p_coupling', float),
        ('p_diffusion', float),
        ('radius_of_gyration', float),
        ('neighbour_degrees', object),  # Keys from the Counter
        ('neighbour_frequencies', object)  # Values from the Counter
    ])

    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            p_c = A[i, j]
            p_d = B[i, j]

            monomer_params = ['A', p_d, 0, 1.0, 0.00, p_c, 0] 
            # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy

            width = 30 # only even numbers

            lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
            # monomers = [Monomer(*monomer_params) for _ in range(50)]

            neighbour_freq, radius, radius_of_gyration = slow_growth_simulation(lattice, monomer_params, total_monomers=20, max_steps=1e6)

            # Extract keys and values from the Counter
            degrees = list(neighbour_freq.keys())
            frequencies = list(neighbour_freq.values())

            results[i, j] = (p_c, p_d, radius_of_gyration, degrees, frequencies)

            print(results[i, j])

    save_data(results, r"C:\Users\User\Desktop\research-updates\thesis\figures\results\2D_KMC\simulation_results.txt")
    print("Results saved to simulation_results.txt.")


"""     fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_surface(A, B, results, cmap = 'viridis')

    ax.set_xlabel("Coupling Probability")
    ax.set_ylabel("Diffusion Probability")
    ax.set_zlabel("Radius of Gyration")

    plt.show() """
    
    # place the monomers at initial positions
    # lattice.randomly_place_monomers(monomers)

    # call the plot_simulation function to visualize the diffusion
    # plot_simulation(lattice, monomers, max_steps = 1000, animate=False)

def save_data(results, filename):
    """
    Save the results matrix to a text file.

    Parameters:
        results (np.ndarray): The structured results array.
        filename (str): The name of the output file.
    """
    with open(filename, 'w') as f:
        # Write a header
        f.write("p_coupling, p_diffusion, radius_of_gyration, neighbour_degrees, neighbour_frequencies\n")

        # Write each entry in the results array
        for i in range(results.shape[0]):
            for j in range(results.shape[1]):
                A_ij = results[i, j]['p_coupling']
                B_ij = results[i, j]['p_diffusion']
                radius = results[i, j]['radius_of_gyration']
                degrees = results[i, j]['neighbour_degrees']
                frequencies = results[i, j]['neighbour_frequencies']
                
                # Format degrees and frequencies as strings
                degrees_str = ','.join(map(str, degrees))
                frequencies_str = ','.join(map(str, frequencies))

                # Write to the file
                f.write(f"{A_ij}, {B_ij}, {radius}, [{degrees_str}], [{frequencies_str}]\n")

def save_meshgrid(A, B, filename_A="meshgrid_p_c.npy", filename_B="meshgrid_p_d.npy"):
    """
    Save meshgrid parts A and B to separate .npy files.

    Parameters:
        A (np.ndarray): The first part of the meshgrid.
        B (np.ndarray): The second part of the meshgrid.
        filename_A (str): Filename for saving A.
        filename_B (str): Filename for saving B.
    """
    np.save(filename_A, A)
    np.save(filename_B, B)


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

def introduce_new_monomer(lattice, new_monomer, monomers, max_steps=1e5):
    '''
    Introduces a new monomer on the lattice and executes the Monomer.action() method iteratively until 
    monomer is coupled (at which point it is defined to not move anymore) or until max_steps has been reached. 
    In the latter case, the monomer is removed from the lattice. This was done to avoid having several islands growing at the same time. 
    We might want to relax this condition at some point.
    '''
    steps = 0
    while steps < max_steps:
        new_monomer.action(lattice)
        if new_monomer.coupled:
            monomers.append(new_monomer)
            print(f"Monomer succesfully coupled after {steps} steps")
            break
        steps += 1
    else:
        """ 
        x, y = new_monomer.get_position()
        lattice.remove_monomer(x, y)
        """
        print(f"Monomer failed to couple after {max_steps} steps. Initializing new monomer...")

def slow_growth_simulation(lattice, monomer_params, total_monomers, max_steps=1e5):
    '''
    What I call slow_growth here is what I had explained in our meeting, where the density of monomers is so low that 
    it is physically accurate to model only a single monomer at a time until it coupled to the growing island. Only after 
    the monomer has coupled is the next one introduced.
    '''
    monomer_1, monomer_2 = initialize_dimer(lattice, monomer_params)
    monomers = [monomer_1, monomer_2]

    for i in range(2, total_monomers):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer]) # initialize monomer with random position (note that this can also be inside the island on an unoccupied site)
        introduce_new_monomer(lattice, new_monomer, monomers, max_steps)
        

    print("Growth simulation completed.")

    neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)

    return neighbour_freq, radius, radius_of_gyration

    # plot_analysis_results(neighbour_freq, radius, lattice, monomers) # some preliminary analysis of the resulting structure

if __name__ == "__main__":
    main()