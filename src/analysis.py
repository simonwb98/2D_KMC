from collections import Counter
import numpy as np
from scipy.spatial.distance import pdist, squareform
import networkx as nx
from skimage.morphology import skeletonize
from skimage.measure import regionprops, label

def analyze_structure(lattice, monomers):
    """Perform analysis on the resulting structure after growth."""
    
    neighbour_frequencies = Counter()
    for monomer in monomers:
        coupled_neighbours = sum(1 for neighbour in lattice.get_next_nearest_neighbours(*monomer.get_position(), monomer.get_orientation()) 
                                if lattice.is_occupied(*neighbour) and lattice.grid[neighbour[1]][neighbour[0]].coupled)
        neighbour_frequencies[coupled_neighbours] += 1
    
    print("Frequency of coupled neighbours:")
    for neighbours, count in sorted(neighbour_frequencies.items()):
        print(f"{count} monomers have {neighbours} coupled neighbours.")

    radius, radius_of_gyration = calculate_effective_radius(lattice, monomers)
    print(f"The Radius of Gyration of the Structure is {radius_of_gyration}")

    # Skeletonize and analyze enclosed areas
    num_enclosed_areas, avg_enclosed_area, enclosed_areas = skeletonize_and_analyze(lattice)
    print(f"Number of enclosed areas: {num_enclosed_areas}")
    print(f"Average enclosed area: {avg_enclosed_area:.2f}")

    # Calculate MST metrics
    positions = np.array([monomer.position for monomer in monomers])
    cell_areas = calculate_average_cell_area(lattice, neighbour_frequencies)
    m, sigma = calculate_mst_metrics(positions, cell_areas)
    print(f"Normalized mean edge length (m): {m:.4f}")
    print(f"Normalized standard deviation of edge lengths (sigma): {sigma:.4f}")

    return neighbour_frequencies, radius, radius_of_gyration, m, sigma, num_enclosed_areas, avg_enclosed_area, enclosed_areas

def skeletonize_and_analyze(lattice):
    """
    Skeletonize the monomer network on the lattice and calculate enclosed area statistics.

    Parameters:
        lattice: Lattice object containing the monomer network.

    Returns:
        num_enclosed_areas (int): Number of enclosed areas.
        average_area (float): Average area of the enclosed regions.
        enclosed_areas (list): List of areas of each enclosed region.
    """
    # Step 1: Convert lattice to binary grid
    binary_grid = np.array([[1 if cell is not None else 0 for cell in row] for row in lattice.grid])

    # Step 2: Skeletonize the binary grid
    skeleton = skeletonize(binary_grid)

    # Step 3: Label enclosed areas
    labeled_regions = label(~skeleton)  # Invert skeleton to find enclosed regions
    props = regionprops(labeled_regions)

    # Step 4: Calculate areas
    enclosed_areas = [region.area for region in props if region.area > 1]  # Exclude trivial areas
    num_enclosed_areas = len(enclosed_areas)
    average_area = np.mean(enclosed_areas) if enclosed_areas else 0

    return num_enclosed_areas, average_area, enclosed_areas

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
    radius_of_gyration = radius
    
    return radius, radius_of_gyration

def calculate_mst_metrics(positions, cell_areas):
    """
    Calculate MST metrics (m, \sigma) from molecular positions and cell areas.

    Parameters:
        positions (ndarray): Array of shape (N, 2) containing the x, y coordinates of molecular centers.
        cell_areas (float): Average cell area in the system.

    Returns:
        m (float): Normalized mean MST edge length.
        sigma (float): Normalized standard deviation of MST edge lengths.
    """
    # Calculate pairwise distances
    distance_matrix = squareform(pdist(positions))

    # Build graph and calculate the minimum spanning tree
    G = nx.from_numpy_array(distance_matrix)
    mst = nx.minimum_spanning_tree(G)

    # Extract MST edge lengths
    edge_lengths = [data['weight'] for _, _, data in mst.edges(data=True)]
    
    # Compute mean and standard deviation of edge lengths
    mean_edge_length = np.mean(edge_lengths)
    std_edge_length = np.std(edge_lengths)

    # Normalize values
    n = len(positions)  # Number of nodes (cells)
    normalization_factor = np.sqrt(cell_areas) * ((n - 1) / n)

    m = mean_edge_length / normalization_factor
    sigma = std_edge_length / normalization_factor

    return m, sigma

def calculate_average_cell_area(lattice, neighbour_frequencies):
    """
    Calculate the average cell area based on lattice geometry and the number of cells.

    Parameters:
        lattice: The lattice object representing the simulation grid.
        neighbour_frequencies: Counter object with frequencies of coupled neighbours.

    Returns:
        float: Average cell area.
    """
    # Number of cells is estimated by the number of monomers with exactly 3 coupled neighbours
    num_cells = neighbour_frequencies.get(3, 0)

    if num_cells == 0:
        raise ValueError("No cells with exactly 3 coupled neighbours found. Cannot calculate cell area.")

    if lattice.rotational_symmetry == 4:
        # Square symmetry: Effective cell area
        effective_cell_area = 1.0  # Assuming unit spacing
    elif lattice.rotational_symmetry == 6:
        # Hexagonal symmetry: Effective cell area
        effective_cell_area = (3 ** 0.5) / 2  # Area of a hexagonal cell with unit spacing
    else:
        raise NotImplementedError("Cell area calculation is only implemented for 4-fold and 6-fold rotational symmetries.")

    return effective_cell_area * num_cells
