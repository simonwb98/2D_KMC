# src/lattice.py
import random
import numpy as np

from monomer import Monomer
from defect import Defect
import math


class Lattice:
    def __init__(self, width, rotational_symmetry, wall=[(0,0), 0, 0], periodic = False, temperature = 300):
        self.width = width
        self.height = width
        self.rotational_symmetry = rotational_symmetry
        self.periodic = periodic
        self.grid = None # will be defined below
        self.lattice_coord = []
        self.temperature = temperature # in K
        self.substrate_properties = {}
        self.define_grid()
        self.wall = wall
        if wall:
            # Wall should be a list of length 3. The first parameter determines the direction of the wall "horiz" or "vert", the second sets the position,
            # the third determines the fraction of the original movement probability for a monomer to move away from a wall
            # the fourth determines the fraction of the original movement probability for a monomer to move over a wall
            # ex) wall=["horiz", 3, 0.3, 0.1]
            self.wall_params = wall
            self.wall_grid = self.make_wall(self.wall_params[0], self.wall_params[1])

    def define_grid(self):
        '''
        Define the grid of the lattice. This does not (yet) include nucleation sites or the spacing between atoms.
        '''
        self.grid, self.lattice_coord = self.generate_grid()
        self.lattice_coord = set(self.lattice_coord)

    def generate_grid(self):
        grid = [[None for i in range(self.width)] for i in range(self.height)]
        
        # Build grid for defects (I made it just create the class parameter here bc adding it to the outputs creates issues elsewhere)
        self.dgrid = [[None for i in range(self.width)] for i in range(self.height)]
        
        lattice_coord = [(i, j) for i in range(self.width) for j in range(self.height)]
        return grid, lattice_coord

    def make_wall(self, direction, pos):
        '''
        Returns a matrix with all zeros, with a 'wall'. 
        The wall is either horizontal or verticle and is made of a line of 2's surrounded by 1's.
        Input the type of wall (either "vert" or "horiz"), and the x or y position of the wall. 
        '''
        m_default = np.zeros((self.width,self.height))

        # check if inputs are in range
        if self.width<pos+1:
            raise ValueError("x0 coordinate not within the size of the grid")
        if pos<0:
            raise ValueError("x0 coordinate not within the size of the grid")
        if self.height<pos+1:
            raise ValueError("y0 coordinate not within the size of the grid")
        if pos<0:
            raise ValueError("y0 coordinate not within the size of the grid")    
        
        if direction == "vert":
            for i in range(self.height):
                m_default[i][pos]=2
            for i in range(self.height):
                for j in range(self.width):
                    if (m_default[i][j]==2):
                        if (pos<self.width-1):
                            m_default[i][j+1]=1

        if direction == "horiz":
            for i in range(self.width):
                m_default[pos][i]=2
            for i in range(self.width):
                for j in range(self.height):
                    if (m_default[i][j]==2):
                        if (pos<self.height-1):
                            m_default[i+1][j]=1

        return m_default

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
    
    # Checks defect grid for defect
    def has_defect(self, x, y):
        return self.dgrid[y][x] is not None
             
    def place_monomer(self, monomer, x, y):
        # first, check if coordinates are properly wrapped
        x, y = self.wrap_coordinates(x, y)

        if not self.is_occupied(x, y):
            self.grid[y][x] = monomer
            monomer.set_position(x, y)
        else:
            pass
    
    def construct_herringbone(self, start_x, start_y, tau, width, nucleation_prob):
        max_length = int(self.width/tau)
        defects = [Defect(0, 0, nucleation_prob)]
        self.place_defect(defects[0], start_x, start_y)
        new_x = start_x
        new_y = start_y
        
        for i in range(max_length*2):
            defects.append(Defect(0, 0, nucleation_prob))
            new_x += math.ceil(tau/2)
            if i % 2 == 0:
                new_y += width - 1
            else:
                new_y -= width - 1
            
            self.place_defect(defects[-1], new_x, new_y)
        
        return defects
        
    
    # Identical to place_monomer, just for defect grid
    def place_defect(self, defect, x, y):
        # first, check if coordinates are properly wrapped
        x, y = self.wrap_coordinates(x, y)

        if not self.has_defect(x, y):
            self.dgrid[y][x] = defect
            defect.set_position(x, y)
        else:
            pass
        
    def remove_defect(self, x, y):
        if self.has_defect(x, y):
            self.dgrid[y][x] = None

    def randomly_place_monomers(self, monomers):
        for monomer in monomers:
            unoccupied = [(x, y) for (x, y) in self.lattice_coord if not self.is_occupied(x, y)]

            if unoccupied:
                (x, y) = random.choice(unoccupied)
                monomer.set_position(x, y)
                self.grid[y][x] = monomer
 
    # Identical to randomly_place_monomers, just for defects
    def randomly_place_defects(self, defects):
        for defect in defects:
            unoccupied = [(x, y) for (x, y) in self.lattice_coord if not self.has_defect(x, y)]

            if unoccupied:
                (x, y) = random.choice(unoccupied)
                defect.set_position(x, y)
                self.dgrid[y][x] = defect
    


    def remove_monomer(self, x, y):
        if self.is_occupied(x, y):
            self.grid[y][x] = None

    def move_monomer(self, monomer, x_new, y_new):
        # moves monomer from old coordinates to new ones
        x_old, y_old = monomer.get_position()
        if not self.is_occupied(x_new, y_new) and not self.has_defect(x_new,y_new):
            self.place_monomer(monomer, x_new, y_new)
            self.remove_monomer(x_old, y_old)
            monomer.set_position(x_new, y_new)
        else:
            pass
    
    # Identical to move_monomer, just for defects
    def move_defect(self, defect, x_new, y_new):
        # moves monomer from old coordinates to new ones
        x_old, y_old = defect.get_position()
        if not self.has_defect(x_new, y_new):
            self.place_defect(defect, x_new, y_new)
            self.remove_defect(x_old, y_old)
            defect.set_position(x_new, y_new)
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

    def get_neighbours_with_wall(self, x, y):
        if self.rotational_symmetry == 4:
            # I added to the output the values of the wall grid for each neighbour point
            neighbours = [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
            n_wrap = list(map(lambda coord: self.wrap_coordinates(*coord), neighbours))
            walls = [self.wall_grid[n_wrap[0][0]][n_wrap[0][1]], self.wall_grid[n_wrap[1][0]][n_wrap[1][1]], self.wall_grid[n_wrap[2][0]][n_wrap[2][1]], self.wall_grid[n_wrap[3][0]][n_wrap[3][1]]]
        elif self.rotational_symmetry == 6:
            if y % 2 == 0:
                neighbours = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y - 1), (x - 1, y + 1)]
                n_wrap = list(map(lambda coord: self.wrap_coordinates(*coord), neighbours))
                walls = [self.wall_grid[n_wrap[0][0]][n_wrap[0][1]], self.wall_grid[n_wrap[1][0]][n_wrap[1][1]], self.wall_grid[n_wrap[2][0]][n_wrap[2][1]], self.wall_grid[n_wrap[3][0]][n_wrap[3][1]], self.wall_grid[n_wrap[4][0]][n_wrap[4][1]], self.wall_grid[n_wrap[5][0]][n_wrap[5][1]]]
            else:
                neighbours = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x + 1, y - 1), (x + 1, y + 1)] 
                n_wrap = list(map(lambda coord: self.wrap_coordinates(*coord), neighbours))
                walls = [self.wall_grid[n_wrap[0][0]][n_wrap[0][1]], self.wall_grid[n_wrap[1][0]][n_wrap[1][1]], self.wall_grid[n_wrap[2][0]][n_wrap[2][1]], self.wall_grid[n_wrap[3][0]][n_wrap[3][1]], self.wall_grid[n_wrap[4][0]][n_wrap[4][1]], self.wall_grid[n_wrap[5][0]][n_wrap[5][1]]]
        return [n_wrap, walls]
        
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
