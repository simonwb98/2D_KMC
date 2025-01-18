# src/monomer.py

import math, random
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
from defect import Defect
class Monomer:
    def __init__(self, monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy, orientations = [0, 180]):
        self.monomer_type = monomer_type
        self.diffusion_rate = diffusion_rate          
        self.diffusion_energy = diffusion_energy
        self.rotation_rate = rotation_rate
        self.rotation_energy = rotation_energy
        self.coupling_rate = coupling_rate
        self.coupling_energy = coupling_energy
        self.coupled = False
        self.nucleating = False
        self.position = None
        self.orientations = orientations
        self.orientation = random.choice(orientations)

    def set_position(self, x, y):
        self.position = (x, y)

    def get_position(self):
        return self.position
    
    def get_orientation(self):
        return self.orientation
    
    def set_orientation(self, orientation):
        self.orientation = orientation
    
    def diffusion_probability(self, lattice):
        '''
        You could modify this to be a property of the lattice class and then implement your matrix model :)
        '''
        temperature = lattice.temperature
        return self.diffusion_rate * math.exp(-self.diffusion_energy / (k_B * temperature)) if not self.coupled else 0 # for two or more coupled monomers, the diffusion probability is (for now) set to 0.
    
    def rotation_probability(self, lattice):
        temperature = lattice.temperature
        return self.rotation_rate * math.exp(-self.rotation_energy / (k_B * temperature)) if not self.coupled else 0
    
    def coupling_probability(self, lattice):
        temperature = lattice.temperature
        return self.coupling_rate * math.exp(-self.coupling_energy / k_B * temperature) if not self.coupled else 0
    
    def couple_with(self, other):
        self.coupled = True
        other.coupled = True


    def diffuse(self, lattice, first_time):
        diffusion_prob = self.diffusion_probability(lattice) # get probability of moving
        # the try block is in case no walls are built
        try:
            grid = lattice.wall_grid
            divisor = lattice.rotational_symmetry
            site0 = lattice.wall_params[2]
            site2 = lattice.wall_params[3]

            x_cur, y_cur = self.get_position()
            if grid[y_cur][x_cur] == 1:
                neighbours_grid = lattice.get_neighbours_with_wall(*self.get_position())[0] # get the coords of neighbours
                neighbours_wall = lattice.get_neighbours_with_wall(*self.get_position())[1] # get the strengths of neighbouring cells
                num_ones = neighbours_wall.count(1)                
                num_twos = neighbours_wall.count(2)
                num_other = divisor - num_twos - num_ones
                index_ones = [i for i, x in enumerate(neighbours_wall) if (x == 1)]
                index_twos = [i for i, x in enumerate(neighbours_wall) if x == 2]
                index_other = [i for i, x in enumerate(neighbours_wall) if (x == 0)]

                p_over_wall = num_twos*site2*diffusion_prob/divisor # probability of jumping over the wall
                p_zero = num_other*site0*diffusion_prob/divisor # probability of moving to a 0
                p_ones = num_ones*diffusion_prob/divisor # probability of moving along the wall
                p_stay = 1 - (p_over_wall + p_zero + p_ones)
                probabilities = [p_zero, p_ones, p_over_wall, p_stay]
                prob_index = random.choices(range(len(probabilities)), weights=probabilities)[0] # choose a probability and return its index

                landing_spots = []
                if prob_index == 0:
                    for index in index_other:
                        landing_spots.append(neighbours_grid[index])
                    x_new, y_new = random.choice(landing_spots)
                    lattice.move_monomer(self, x_new, y_new)

                elif prob_index == 1:
                    for index in index_ones:
                        landing_spots.append(neighbours_grid[index])
                    x_new, y_new = random.choice(landing_spots)
                    lattice.move_monomer(self, x_new, y_new)
                    
                elif prob_index == 2:
                    for index in index_twos:
                        landing_spots.append(neighbours_grid[index])
                    x_new, y_new = random.choice(landing_spots)
                    lattice.move_monomer(self, x_new, y_new)

                else:
                    if first_time:
                        self.coupled = True
                    else:
                        self.couple_with_wall(lattice)

            if grid[y_cur][x_cur] == 2:
                neighbours_grid = lattice.get_neighbours_with_wall(*self.get_position())[0] # get the coords of neighbours
                neighbours_wall = lattice.get_neighbours_with_wall(*self.get_position())[1] # get the strengths of neighbouring cells
                num_ones = neighbours_wall.count(1)
                num_twos = neighbours_wall.count(2)
                num_other = divisor - num_twos - num_ones
                index_ones = [i for i, x in enumerate(neighbours_wall) if (x == 1)]
                index_twos = [i for i, x in enumerate(neighbours_wall) if x == 2]
                index_other = [i for i, x in enumerate(neighbours_wall) if (x == 0)]

                p_over_wall = num_ones*site0*diffusion_prob/divisor # probability of jumping over the wall
                p_zero = num_other*site0*diffusion_prob/divisor # probability of moving to a 0
                p_twos = num_twos*diffusion_prob/divisor # probability of moving along the wall
                p_stay = 1 - (p_over_wall + p_zero + p_twos)
                probabilities = [p_zero, p_over_wall, p_twos, p_stay]
                prob_index = random.choices(range(len(probabilities)), weights=probabilities)[0] # choose a probability and return its index

                landing_spots = []
                if prob_index == 0:
                    for index in index_other:
                        landing_spots.append(neighbours_grid[index])
                    x_new, y_new = random.choice(landing_spots)
                    lattice.move_monomer(self, x_new, y_new)

                elif prob_index == 1:
                    for index in index_ones:
                        landing_spots.append(neighbours_grid[index])
                    x_new, y_new = random.choice(landing_spots)
                    lattice.move_monomer(self, x_new, y_new)

                elif prob_index == 2:
                    for index in index_twos:
                        landing_spots.append(neighbours_grid[index])
                    x_new, y_new = random.choice(landing_spots)
                    lattice.move_monomer(self, x_new, y_new)

            else:
                if random.random() < diffusion_prob and not self.nucleating and not self.coupled: # based on the probability, decide if diffuse or not
                    neighbours = lattice.get_neighbours(*self.get_position())
                    x_new, y_new = random.choice(neighbours)
                    lattice.move_monomer(self, x_new, y_new)        
            
     
        except AttributeError:
            if random.random() < diffusion_prob and not self.nucleating and not self.coupled: # based on the probability, decide if diffuse or not
                neighbours = lattice.get_neighbours(*self.get_position())
                x_new, y_new = random.choice(neighbours)
                lattice.move_monomer(self, x_new, y_new)
            

    def rotate(self, lattice):
        rotation_prob = self.rotation_probability(lattice)

        if random.random() < rotation_prob:
            self.set_orientation(random.choice([o for o in self.orientations if not o == self.orientation]))

    def couple(self, lattice):
        coupling_prob = self.coupling_probability(lattice)

        if random.random() < coupling_prob:
            neighbours = lattice.get_neighbours(*self.get_position())
            if any(lattice.is_occupied(nx, ny) for (nx, ny) in neighbours) or lattice.has_defect(*self.get_position()):
                return # disallow coupling when there are nearest neighbours to this monomer. This condition basically realizes the fact that monomers physically restrict each other (geometric hindrance) - a cleaner way of doing this is to disallow diffusion into sites that have monomers that are nearest neighbours, but this is fine also.
            next_neighbours = lattice.get_next_nearest_neighbours(*self.get_position(), self.get_orientation()) # this could potentially be made faster sometime down the line
            neighbouring_monomers = [lattice.grid[ny][nx] for (nx, ny) in next_neighbours if not lattice.grid[ny][nx] == None and ((lattice.grid[ny][nx].get_orientation() != self.get_orientation()) and not isinstance(lattice.dgrid[ny][nx], Defect))]
            
            
            neighboring_defects = [lattice.dgrid[ny][nx] for (nx, ny) in next_neighbours if not lattice.dgrid[ny][nx] == None and isinstance(lattice.dgrid[ny][nx], Defect)]
            if neighbouring_monomers:
                partner = random.choice(neighbouring_monomers)
                partner_neighbours = lattice.get_neighbours(*partner.get_position())
                self.couple_with(partner)
            elif neighboring_defects:
                choice = random.choice(neighboring_defects)
                x, y = choice.get_position()
                
                # Checks if the defect already has a monomer on it
                if choice.nucleating:
                    return
                
                # Move the monomer on top of the defect
                new_x, new_y = choice.get_position()
                self.set_position(new_x, new_y)
                
                self.nucleating = True
                choice.nucleating = True
                            
                return 

    def couple_with_wall(self, lattice):
        coupling_prob = 1

        if random.random() < coupling_prob:
            neighbours = lattice.get_neighbours(*self.get_position())
            if any(lattice.is_occupied(nx, ny) for (nx, ny) in neighbours):
                return # disallow coupling when there are nearest neighbours to this monomer. This condition basically realizes the fact that monomers physically restrict each other (geometric hindrance) - a cleaner way of doing this is to disallow diffusion into sites that have monomers that are nearest neighbours, but this is fine also.
            next_neighbours = lattice.get_next_nearest_neighbours(*self.get_position(), self.get_orientation()) # this could potentially be made faster sometime down the line
            neighbouring_monomers = [lattice.grid[ny][nx] for (nx, ny) in next_neighbours if not lattice.grid[ny][nx] == None and lattice.grid[ny][nx].get_orientation() != self.get_orientation()]

            if neighbouring_monomers:
                partner = random.choice(neighbouring_monomers)
                partner_neighbours = lattice.get_neighbours(*partner.get_position())

                if any(lattice.is_occupied(nx, ny) for (nx, ny) in partner_neighbours):
                    return

                self.coupled = True
                     
    def action(self, lattice, first_time):
        '''
        Runs all actions a monomer can perform in succession.
        '''
        self.diffuse(lattice, first_time)
        self.rotate(lattice)
        self.couple(lattice)
