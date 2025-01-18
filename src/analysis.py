# src/analysis.py

from collections import Counter
import numpy as np

def analyze_structure(lattice, monomers):
    """Perform analysis on the resulting structure after growth."""
    
    neighbour_frequencies = Counter()
    for monomer in monomers:
        coupled_neighbours = sum(1 for neighbour in lattice.get_next_nearest_neighbours(*monomer.get_position(), monomer.get_orientation()) 
                                if lattice.is_occupied(*neighbour) and (lattice.grid[neighbour[1]][neighbour[0]].coupled))
        neighbour_frequencies[coupled_neighbours] += 1
    
    print("Frequency of coupled neighbours:")
    for neighbours, count in sorted(neighbour_frequencies.items()):
        print(f"{count} monomers have {neighbours} coupled neighbours.")

    radius, radius_of_gyration = calculate_effective_radius(lattice, monomers)
    print(f"The Radius of Gyration of the Structure is {radius_of_gyration}")
    
    return neighbour_frequencies, radius, radius_of_gyration

def calculate_effective_radius(lattice, monomers):
    """Calculate the effective radius (radius of gyration) of the resulting structure."""
    
    # calculate the center of mass
    x_positions = np.array([monomer.position[0] for monomer in monomers])
    y_positions = np.array([monomer.position[1] for monomer in monomers])
    
    x_center_of_mass = np.mean(x_positions)
    y_center_of_mass = np.mean(y_positions)
    
    # calculate the squared distances from the center of mass
    squared_distances = (x_positions - x_center_of_mass)**2 + (y_positions - y_center_of_mass)**2
    
    # compute the radius of gyration (effective radius)
    radius = np.sqrt(np.mean(squared_distances))
    radius_of_gyration = radius/np.sqrt(len(monomers))
    
    return radius, radius_of_gyration
