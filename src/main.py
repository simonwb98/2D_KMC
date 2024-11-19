# src/main.py

from lattice import Lattice
from monomer import Monomer
from plotter import plot_simulation, plot_final_state, plot_analysis_results
from analysis import analyze_structure
from kinetic_monte_carlo import kmc_simulation
import random
import numpy as np
from matplotlib import pyplot as plt

### KMC

def initialize_dimer(lattice, monomer_params):
    """
    Initialize a dimer in the center of the lattice.

    Args:
        lattice (Lattice): The lattice object.
        monomer_params (dict): Parameters for monomers.

    Returns:
        list: A list containing the two dimer monomers.
    """
    # Determine center of the lattice
    x_center = lattice.width // 2
    y_center = lattice.height // 2

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

    return [monomer_1, monomer_2]

def grow_dimer(lattice, monomer_params, total_monomers, max_steps=1e5):
    """
    Grow a dimer by sequentially adding monomers until they are coupled.

    Args:
        lattice (Lattice): The lattice object.
        monomer_params (dict): Parameters for monomers.
        total_monomers (int): Total number of monomers to add.
        max_steps (int): Maximum number of kMC steps.

    Returns:
        float: Total simulation time.
    """
    # Initialize the dimer
    dimer = initialize_dimer(lattice, monomer_params)

    # Add monomers sequentially
    monomers = dimer
    total_time = 0

    for _ in range(total_monomers - 2):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer])
        monomers.append(new_monomer)

        # Run kMC simulation for this monomer to couple
        time_spent = kmc_simulation(lattice, monomers, max_steps=max_steps)
        total_time += time_spent

        # Check if the monomer coupled; if not, remove it
        if not new_monomer.coupled:
            lattice.remove_monomer(*new_monomer.get_position())
            monomers.remove(new_monomer)

    print("Growth simulation completed.")
    return lattice, monomers

def main():
    """
    Main method to run the kinetic Monte Carlo simulation.
    """
    # Initialize the lattice
    lattice = Lattice(width=20, rotational_symmetry=6, periodic=True)

    # Create monomers with example parameters
    monomer_params = ["A", 1e12, 0.3, 1e12, 1, 1e12, 0.5]
    # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy


    # Optional: Grow dimer for demonstration
    total_monomers = 30
    max_steps = 100_000
    print("Starting dimer growth simulation...")
    lattice, monomers = grow_dimer(lattice, monomer_params, total_monomers=total_monomers, max_steps=max_steps)
    plot_final_state(lattice, monomers)
    # print(f"Dimer growth simulation completed in {total_time:.2f} units.")

if __name__ == "__main__":
    main()
""" 

""" ### Nucleation Sites
""" 
def main():
    p_c = 0.1
    p_d = 1

    monomer_params = ['A', p_d, 0, 1.0, 0.00, p_c, 0] 
    # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy

    width = 80 # only even numbers

    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
    # monomers = [Monomer(*monomer_params) for _ in range(50)]

    neighbour_freq, radius, radius_of_gyration = nucleation_site_simulation(lattice, monomer_params, total_monomers=100, max_steps=1e6, nucleation_site_density=1/(10 * width))
    
     """


### Copolymers

""" def main():
    p_c_A, p_c_B = (1, 1)
    p_d_A, p_d_B = (1, 1)

    monomer_params_A = ['A', p_d_A, 0, 1.0, 0.00, p_c_A, 0] 
    monomer_params_B = ['B', p_d_B, 0, 1.0, 0.00, p_c_B, 0]
    # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy

    width = 30 # only even numbers

    lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
    # monomers = [Monomer(*monomer_params) for _ in range(50)]

    neighbour_freq, radius, radius_of_gyration = slow_growth_simulation(lattice, monomer_params, total_monomers=20, max_steps=1e6)
 """
    
### Batch Analysis for Different Coupling and Diffusion Probabilities


""" def main():
    N = 2 # N^2 is the number of simulations

    p_c = np.linspace(0.1, 1, N)
    p_d = np.linspace(0.1, 1, N)

    A, B = np.meshgrid(p_c, p_d)
    save_meshgrid(A, B, 
                  r"data\meshgrid_p_c.npy", 
                  r"data\meshgrid_p_d.npy")
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

    save_data(results, r"data\simulation_results.txt")
    print("Results saved to simulation_results.txt.") """


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

""" def save_data(results, filename):
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
                f.write(f"{A_ij}, {B_ij}, {radius}, [{degrees_str}], [{frequencies_str}]\n") """

""" def save_meshgrid(A, B, filename_A="meshgrid_p_c.npy", filename_B="meshgrid_p_d.npy"):
    
    np.save(filename_A, A)
    np.save(filename_B, B) """


# Approaching a more phyisically meaningful system - Reducing the complexity by adding one monomer at a time

""" def initialize_dimer(lattice, monomer_params, x=0, y=0):
    x_center, y_center = (x, y) # find centre
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
 """
""" def introduce_new_monomer(lattice, new_monomer, monomers, max_steps=1e5):
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
        
        print(f"Monomer failed to couple after {max_steps} steps. Initializing new monomer...") """

""" def slow_growth_simulation(lattice, monomer_params, total_monomers, max_steps=1e5):
    monomer_1, monomer_2 = initialize_dimer(lattice, monomer_params, x=lattice.width//2, y=lattice.width//2)
    monomers = [monomer_1, monomer_2]

    for i in range(2, total_monomers):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer]) # initialize monomer with random position (note that this can also be inside the island on an unoccupied site)
        introduce_new_monomer(lattice, new_monomer, monomers, max_steps)
        

    print("Growth simulation completed.")

    neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)

    return neighbour_freq, radius, radius_of_gyration
 """
    # plot_analysis_results(neighbour_freq, radius, lattice, monomers) # some preliminary analysis of the resulting structure

""" def calculate_hexagonal_nucleation_sites(lattice, nucleation_site_density, boundary_offset=2):
    num_sites = int(nucleation_site_density * (lattice.width - 2 * boundary_offset) ** 2)
    spacing = int(((lattice.width - 2 * boundary_offset) * (lattice.height - 2 * boundary_offset) / num_sites) ** 0.5)

    coordinates = []
    for j in range(boundary_offset, lattice.height - boundary_offset, spacing):  # y spacing with sqrt(3)/2 factor
        for i in range(boundary_offset, lattice.width - boundary_offset, spacing):
            x = i + (spacing // 2 if (j // spacing) % 2 else 0)  # Offset every other row
            if x < lattice.width - boundary_offset and j < lattice.height - boundary_offset:
                coordinates.append((x, j))

    return coordinates """


""" def nucleation_site_simulation(
    lattice, monomer_params, total_monomers, max_steps=1e5, nucleation_site_density=1/300
):
    
    # Step 1: Calculate nucleation site coordinates in hexagonal arrangement
    nucleation_sites = calculate_hexagonal_nucleation_sites(lattice, nucleation_site_density, boundary_offset=lattice.width//10)
    print(f"Calculated {len(nucleation_sites)} nucleation sites in a hexagonal arrangement.")

    # Step 2: Initialize dimers at the calculated positions
    dimers = []
    for x, y in nucleation_sites:
        monomer_1, monomer_2 = initialize_dimer(lattice, monomer_params, x=x, y=y)
        dimers.extend([monomer_1, monomer_2])

    # Step 3: Introduce additional monomers for slow growth
    monomers = dimers[:]
    for i in range(len(dimers), total_monomers):
        new_monomer = Monomer(*monomer_params)
        lattice.randomly_place_monomers([new_monomer])  # Random initial position
        introduce_new_monomer(lattice, new_monomer, monomers, max_steps)  # Growth step
        monomers.append(new_monomer)

    print("Growth simulation completed.")

    # Step 4: Analyze the resulting structure
    neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)
    plot_analysis_results(neighbour_freq, radius_of_gyration, lattice, monomers)

    return neighbour_freq, radius, radius_of_gyration


if __name__ == "__main__":
    main()  """