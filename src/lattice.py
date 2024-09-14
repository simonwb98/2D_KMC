# src/lattice.py
import random

class Lattice:
    def __init__(self, width, rotational_symmetry, periodic = False, temperature = 300):
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
        if self.rotational_symmetry == 4:
            self.grid, self.lattice_coord = self.generate_square_grid()
        elif self.rotational_symmetry == 6:
            self.grid, self.lattice_coord = self.generate_hex_grid()
        else:
            raise ValueError(f"Rotational symmetry {self.rotational_symmetry} not supported. Please check input.")
        
        self.lattice_coord = set(self.lattice_coord)

    def generate_square_grid(self):
        grid = [[None for i in range(self.width)] for i in range(self.height)]
        lattice_coord = [(i, j) for i in range(self.width) for j in range(self.height)]
        return grid, lattice_coord
    
    def generate_hex_grid(self):
        grid = []
        lattice_coord = []
        for j in range(self.height):
            if j % 2 == 0: # offset every second row
               row = [None]*self.width
               row_coord = [(i, j) for i in range(self.width)]
            else:
                row = [None]*(self.width - 1)
                row_coord = [(i, j) for i in range(self.width - 1)]
            grid.append(row)
            lattice_coord.extend(row_coord)
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
                        x = self.width - 2
                    elif x >= self.width - 1:  # right to left
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