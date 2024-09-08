# src/polymer_simulation/lattice.py

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
        for i in range(self.height):
            if i % 2 == 0: # offset every second row
               row = [None]*self.width
               row_coord = [(i, j) for j in range(self.width)]
            else:
                row = [None]*(self.width - 1)
                row_coord = [(i, j) for j in range(self.width - 1)]
            grid.append(row)
            lattice_coord.extend(row_coord)
        return grid, lattice_coord
    
    def is_member(self, x, y):
        return (x, y) in self.lattice_coord

    def is_occupied(self, x, y):
        return self.grid[y][x] is not None

    def place_monomer(self, monomer, x, y):
        if not self.is_occupied(x, y):
            self.grid[y][x] = monomer

    def remove_monomer(self, x, y):
        if self.is_occupied(x, y):
            self.grid[y][x] = None

    def get_neighbors(self, x, y):
        if self.rotational_symmetry == 4:
            neighbours = [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
        elif self.rotational_symmetry == 6:
            if y % 2 == 0:
                neighbours = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y - 1), (x - 1, y + 1)]
            else:
                neighbours = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y + 1), (x + 1, y + 1)]

        # periodic boundary conditions
        if self.periodic:
            adjusted_neighbors = []
            for nx, ny in neighbours:
                if not self.is_member(nx, ny):
                    if ny < 0:  # wrap coordinates vertically (top to bottom)
                        ny = self.height - 1
                    elif ny >= self.height:  # bottom to top
                        ny = 0
                    
                    if ny % 2 == 0:
                        if nx < 0:  # wrap horizontally on even rows (left to right)
                            nx = self.width - 1
                        elif nx >= self.width:  # right to left
                            nx = 0
                    else:  
                        if nx < 0:  # wrap horizontally on odd rows (left to right)
                            nx = self.width - 2
                        elif nx >= self.width - 1:  # right to left
                            nx = 0
                
                adjusted_neighbors.append((nx, ny))
            neighbours = adjusted_neighbors
        else:
            neighbours = [(nx, ny) for nx, ny in neighbours if self.is_member(nx, ny)]
        return neighbours
        
