# src/main.py

from lattice import Lattice
from monomer import Monomer
from defect import Defect
from plotter import plot_simulation, plot_final_state, plot_analysis_results
import matplotlib.pyplot as plt
# from plotter import plot_simulation, plot_final_state, plot_analysis_results

from analysis import analyze_structure
import random
import numpy as np

# from matplotlib import pyplot as plt
import csv

# COMMENT THE NEXT FIVE FUNCTIONS OUT IF YOU'D LIKE, JUST DOING THIS FOR tESTING
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

def introduce_new_monomer(lattice, new_monomer, monomers, first_time, defects, max_steps=1e6):

    '''
    Introduces a new monomer on the lattice and executes the Monomer.action() method iteratively until 
    monomer is coupled (at which point it is defined to not move anymore) or until max_steps has been reached. 
    In the latter case, the monomer is removed from the lattice. This was done to avoid having several islands growing at the same time. 
    We might want to relax this condition at some point.
    '''
    steps = 0
    while steps < max_steps:
        new_monomer.action(lattice, first_time)

        for mon in monomers:
            mon.dehalogenate(lattice)
        if new_monomer.coupled or new_monomer.nucleating:
            monomers.append(new_monomer)
            print(f"Monomer succesfully coupled after {steps} steps")
            return 0
        
        steps += 1

    else:
        x, y = new_monomer.get_position()
        lattice.remove_monomer(x, y)
        print(f"Monomer failed to couple after {max_steps} steps. Initializing new monomer...")
        return 1

def create_defects(defect_density, lattice, defect_params):
    num_defects = round(defect_density*lattice.width**2)
    defects = []
    for num in range(num_defects):
        defects.append(Defect(*defect_params))
    return defects

def slow_growth_simulation(lattice, monomer_params, defect_params, defect_density, total_monomers, herringbone=True, max_steps=1e6):
    '''
    What I call slow_growth here is what I had explained in our meeting, where the density of monomers is so low that 
    it is physically accurate to model only a single monomer at a time until it coupled to the growing island. Only after 
    the monomer has coupled is the next one introduced.
    '''

    # Change the initialization of the dimer to a normal introduction of one monomer and allow it to nucleate at some point
    monomer_1, monomer_2 = initialize_dimer(lattice, monomer_params)
    #monomers = [monomer_1, monomer_2]
    monomers = [monomer_1, monomer_2]
    first_time = True
    for i in range(2, total_monomers):
        new_monomer = Monomer(*monomer_params) # Constructing a monomer here is fine, but
                                               # we might want to input the lattice into the monomer as a matrix of probabilities
        j = 1                     
        lattice.randomly_place_monomers([new_monomer]) # initialize monomer with random position (note that this can also be inside the island on an unoccupied site)
        while j==1:
            j = introduce_new_monomer(lattice, new_monomer, monomers, first_time, max_steps)
        
        first_time = False
        

    print("Growth simulation completed.")

    #neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)

    #plot_analysis_results(neighbour_freq, radius, lattice, monomers) # some preliminary analysis of the resulting structure
    return monomers

def save_results_to_csv(results, filename):
    """
    Save aggregated results to a CSV file.

    Args:
        results (list of dict): Aggregated results to save.
        filename (str): Output file name.
    """
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Write header
        header = [
            "Diffusion Energy", 
            "Rotation Energy", 
            "Coupling Energy", 
            "Dehalogenation Energy",
            "Averaged Neighbour Frequency", 
            "Average Radius", 
            "Standard Dev Radius",
            "Average Radius of Gyration",
            "Standard Dev ROG"
        ]
        writer.writerow(header)

        # Write rows
        for result in results:
            writer.writerow([
                result["diffusion_energy"],
                result["rotation_energy"],
                result["coupling_energy"],
                result["dehalogen_energy"],
                result["averaged_neighbour_freq"],
                result["avg_radius"],
                result["std_radius"],
                result["avg_radius_of_gyration"],
                result["std_radius_of_gyration"]
            ])


def main():
    # Initialize lattice and monomers
    diffusion_energies = np.linspace(0,1.5,6)
    coupling_energies = np.linspace(0,1.5,6)
    dehalogen_energies = np.linspace(0, 2.3, 6)
    rotation_energies = [0]

    width = 60 # only even numbers
    num_simulations_per_triplet = 3
    aggregated_results = [] # Store aggregated results
    axes = []
    for diff_energy in diffusion_energies:
        for rot_energy in rotation_energies:
            for coup_energy in coupling_energies:
                for dehal_energy in dehalogen_energies:
                    monomer_params = ['A', 1e13, diff_energy, 1e13, rot_energy, 1e13, coup_energy, 1e13, dehal_energy] 
                    # monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy

                    defect_params = [1.0, 0.00, 1.0]
                    # diffusion_rate, diffusion_energy, nucleation_prob

                    # Initialize containers for results
                    all_neighbour_freqs = []
                    all_radii = []
                    all_radii_of_gyration = []
                    
                    print('###########################################################################################################')
                    print('')
                    print(f"Simulation series beginning with diffusion: {diff_energy}, rotation: {rot_energy}, coupling: {coup_energy}, dehalogen: {dehal_energy}")
                    print('')
                    print('###########################################################################################################')
                    
                    for sim in range(num_simulations_per_triplet):
                        lattice = Lattice(width=width, rotational_symmetry=6, periodic=True)
                        monomers = slow_growth_simulation(lattice, monomer_params,defect_params,defect_density=0.0, total_monomers=50, max_steps=1e6)
        
                        # call the plot_simulation function to visualize the diffusion
                        #plot_simulation(lattice, monomers, max_steps = 1000, animate=False)
                        neighbour_freq, radius, radius_of_gyration = analyze_structure(lattice, monomers)
                        all_neighbour_freqs.append(neighbour_freq)
                        all_radii.append(radius)
                        all_radii_of_gyration.append(radius_of_gyration)
                    #axes.append(plot_analysis_results(neighbour_freq, radius, lattice, monomers))
                    combined_freqs = {}
                    for freq in all_neighbour_freqs:
                        for degree, count in freq.items():
                            combined_freqs[degree] = combined_freqs.get(degree, 0) + count

                    averaged_neighbour_freq = {degree: count / num_simulations_per_triplet for degree, count in combined_freqs.items()}

                    # Average radius and radius of gyration
                    avg_radius = np.mean(all_radii)
                    std_radius = np.std(all_radii)
                    avg_radius_of_gyration = np.mean(all_radii_of_gyration)
                    std_radius_of_gyration = np.std(all_radii_of_gyration)

                    aggregated_results.append({
                        "diffusion_energy": diff_energy,
                        "rotation_energy": rot_energy,
                        "coupling_energy": coup_energy,
                        "dehalogen_energy": dehal_energy,
                        "averaged_neighbour_freq": averaged_neighbour_freq,
                        "avg_radius": avg_radius,
                        "std_radius": std_radius,
                        "avg_radius_of_gyration": avg_radius_of_gyration,
                        "std_radius_of_gyration":std_radius_of_gyration
                    })
                    

    # Save aggregated results to a CSV file
    save_results_to_csv(aggregated_results, r"zach_output_rot.csv")
    
    print("Parameter sweep completed. Results saved to 'zach_output_rot.csv'.")



if __name__ == "__main__":
    a = main()

