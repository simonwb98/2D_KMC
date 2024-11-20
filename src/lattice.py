# src/lattice.py
import random
import numpy as np

class Lattice:
    def __init__(self, width, rotational_symmetry, periodic = False, temperature = 600):
        self.width = width
        self.height = width
        self.rotational_symmetry = rotational_symmetry
        self.periodic = periodic
        self.grid = None # will be defined below
        self.lattice_coord = []
        self.temperature = temperature # in K
        self.substrate_properties = {}

        self.define_grid()

    def define_grid(self):
        '''
        Define the grid of the lattice. This does not (yet) include nucleation sites or the spacing between atoms.
        '''
        self.grid, self.lattice_coord = self.generate_grid()
        self.lattice_coord = set(self.lattice_coord)

    def generate_grid(self):
        grid = [[None for i in range(self.width)] for i in range(self.height)]
        lattice_coord = [(i, j) for i in range(self.width) for j in range(self.height)]
        return grid, lattice_coord
    
    def is_member(self, x, y):
        if (x, y) not in self.lattice_coord:
            raise KeyError(f"Coordinates ({x}, {y}) not found in lattice with lattice sites: {self.lattice_coord}\n")
        return True
        
    def wrap_coordinates(self, x, y):
        '''
        This might be moved into an auxiliary file of helper functions later on. 

        This function wraps the coordinates of an input set of coordinates (x, y) and wraps them around if periodic bc are applied. 
        It takes care of hex and square symmetry grids so far. It has been tested for odd number Lattice.height/Lattice.width.
        '''
        # periodic boundary conditions
        if self.periodic and self.rotational_symmetry == 6:
                if y < 0:  # wrap coordinates vertically (top to bottom)
                    y = self.height - 1
                elif y >= self.height:  # bottom to top
                    y = 0

                if y % 2 == 0:
                    if x < 0:  # wrap horizontally on even rows (left to right)
                        x = self.width - 1
                    elif x >= self.width:  # right to left
                        x = 0
                else:  
                    if x < 0:  # wrap horizontally on odd rows (left to right)
                        x = self.width - 1
                    elif x >= self.width:  # right to left
                        x = 0
        elif self.periodic and self.rotational_symmetry == 4:
            if y < 0:
                y = self.height - 1
            elif y >= self.height:
                y = 0
            
            if x < 0:
                x = self.width - 1
            elif x >= self.width:
                x = 0

        if self.is_member(x, y): # check if coordinates are part of the lattice.
            return (x, y) 
    
    def is_occupied(self, x, y):
        return self.grid[y][x] is not None
             
    def place_monomer(self, monomer, x, y):
        # first, check if coordinates are properly wrapped
        x, y = self.wrap_coordinates(x, y)

        if not self.is_occupied(x, y):
            self.grid[y][x] = monomer
            monomer.set_position(x, y)
        else:
            pass

    def randomly_place_monomers(self, monomers):
        for monomer in monomers:
            unoccupied = [(x, y) for (x, y) in self.lattice_coord if not self.is_occupied(x, y)]

            if unoccupied:
                (x, y) = random.choice(unoccupied)
                monomer.set_position(x, y)
                self.grid[y][x] = monomer
            
    def randomly_place_monomers_at_edge(self, monomers):
        """
        Place monomers randomly at the edge of the lattice.

        The edge is defined as lattice sites along the boundary:
        - Top row (y = 0)
        - Bottom row (y = height - 1)
        - Left column (x = 0)
        - Right column (x = width - 1)

        Args:
            monomers (list): List of monomer objects to place.
        """
        # Define edge coordinates
        edge_coords = []

        # Top and bottom rows
        for x in range(self.width):
            if not self.is_occupied(x, 0):
                edge_coords.append((x, 0))  # Top row
            if not self.is_occupied(x, self.height - 1):
                edge_coords.append((x, self.height - 1))  # Bottom row

        # Left and right columns
        for y in range(1, self.height - 1):  # Avoid duplicates at corners
            if not self.is_occupied(0, y):
                edge_coords.append((0, y))  # Left column
            if not self.is_occupied(self.width - 1, y):
                edge_coords.append((self.width - 1, y))  # Right column

        # Ensure edge sites are unique
        edge_coords = list(set(edge_coords))

        # Place monomers at random edge positions
        for monomer in monomers:
            if edge_coords:
                x, y = random.choice(edge_coords)
                self.place_monomer(monomer, x, y)
                edge_coords.remove((x, y))  # Prevent reuse of the same site
            else:
                raise ValueError("Not enough edge sites to place all monomers.") 

    def remove_monomer(self, x, y):
        if self.is_occupied(x, y):
            self.grid[y][x] = None

    def move_monomer(self, monomer, x_new, y_new):
        # moves monomer from old coordinates to new ones
        x_old, y_old = monomer.get_position()
        if not self.is_occupied(x_new, y_new):
            self.place_monomer(monomer, x_new, y_new)
            self.remove_monomer(x_old, y_old)
            monomer.set_position(x_new, y_new)
        else:
            pass

    def get_neighbours(self, x, y):
        if self.rotational_symmetry == 4:
            neighbours = [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
        elif self.rotational_symmetry == 6:
            if y % 2 == 0:
                neighbours = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y - 1), (x - 1, y + 1)]
            else:
                neighbours = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x + 1, y - 1), (x + 1, y + 1)]

        neighbours = list(map(lambda coord: self.wrap_coordinates(*coord), neighbours))
        return neighbours

        
    def get_next_nearest_neighbours(self, x, y, orientation):
        # for now only 6-fold rotational symmetry
        if self.rotational_symmetry == 6:
            if orientation == 180:
                if y % 2 == 0:
                    next_nearest = [(x, y - 2), (x + 1, y + 1), (x - 2, y + 1)]
                else:
                    next_nearest = [(x, y - 2), (x + 2, y + 1), (x - 1, y + 1)]
            elif orientation == 0:
                if y % 2 == 0:
                    next_nearest = [(x - 2, y - 1), (x + 1, y - 1), (x, y + 2)]
                else:
                    next_nearest = [(x - 1, y - 1), (x + 2, y - 1), (x, y + 2)]
        else:
            raise NotImplementedError("Next nearest neighbour is restricted to 6-fold rotational symmetries (for now).")
        
        next_nearest = list(map(lambda coord: self.wrap_coordinates(*coord), next_nearest))
        return next_nearest
    
    def find_cells(self):
        """
        Identify all enclosed areas (cells) in the lattice.
        
        Returns:
            List of sets: Each set contains the coordinates of monomers that form a cell.
        """
        visited = set()
        cells = []

        def is_cell_boundary(coord):
            """Check if a coordinate is part of a potential cell boundary."""
            x, y = coord
            return self.is_occupied(x, y)

        def explore_boundary(start):
            """Explore the boundary of a cell starting from a given coordinate."""
            stack = [start]
            boundary = set()
            
            while stack:
                current = stack.pop()
                if current not in visited and is_cell_boundary(current):
                    visited.add(current)
                    boundary.add(current)
                    neighbours = self.get_neighbours(*current)
                    stack.extend(neighbours)
            
            return boundary

        for coord in self.lattice_coord:
            if coord not in visited and is_cell_boundary(coord):
                boundary = explore_boundary(coord)
                if boundary:  # A valid enclosed area
                    cells.append(boundary)

        return cells

    
    def precompute_neighbors(self):
        """
        Precompute neighbors and next-nearest neighbors for each lattice site and orientation.
        """
        self.neighbors = {}
        self.next_nearest_neighbors = {}
        orientations = [0, 180]  # Example orientations; adjust based on your system

        for x in range(self.width):
            for y in range(self.height):
                for orientation in orientations:
                    # Compute neighbors and next-nearest neighbors for each orientation
                    self.neighbors[(x, y, orientation)] = self.get_neighbours(x, y)
                    self.next_nearest_neighbors[(x, y, orientation)] = self.get_next_nearest_neighbours(x, y, orientation)

    def get_cached_neighbors(self, x, y, orientation):
        """
        Retrieve precomputed neighbors for a site and orientation.
        """
        return self.neighbors.get((x, y, orientation), [])

    def get_cached_next_nearest_neighbors(self, x, y, orientation):
        """
        Retrieve precomputed next-nearest neighbors for a site and orientation.
        """
        return self.next_nearest_neighbors.get((x, y, orientation), [])


