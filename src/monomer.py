# src/monomer.py

import math, random
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

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

    def diffuse(self, lattice):
        diffusion_prob = self.diffusion_probability(lattice) # get probability
        
        if random.random() < diffusion_prob: # based on the probability, decide if diffuse or not
            neighbours = lattice.get_neighbours(*self.get_position())
            x_new, y_new = random.choice(neighbours)
            lattice.move_monomer(self, x_new, y_new)

    def rotate(self, lattice):
        rotation_prob = self.rotation_probability(lattice)

        if random.random() < rotation_prob:
            self.set_orientation(random.choice([o for o in self.orientations if not o == self.orientation]))

    def couple(self, lattice):
        coupling_prob = self.coupling_probability(lattice)
        neighbours = lattice.get_neighbours(*self.get_position())
        if any(lattice.is_occupied(nx, ny) for (nx, ny) in neighbours):
            return # disallow coupling when there are nearest neighbours to this monomer. This condition basically realizes the fact that monomers physically restrict each other (geometric hindrance) - a cleaner way of doing this is to disallow diffusion into sites that have monomers that are nearest neighbours, but this is fine also.
        next_neighbours = lattice.get_next_nearest_neighbours(*self.get_position(), self.get_orientation()) # this could potentially be made faster sometime down the line
        neighbouring_monomers = [lattice.grid[ny][nx] for (nx, ny) in next_neighbours if not lattice.grid[ny][nx] == None and lattice.grid[ny][nx].get_orientation() != self.get_orientation()]

        if neighbouring_monomers:
            for partner in neighbouring_monomers:
                if random.random() < coupling_prob:
                    partner_neighbours = lattice.get_neighbours(*partner.get_position())

                    if any(lattice.is_occupied(nx, ny) for (nx, ny) in partner_neighbours):
                            return

                    self.couple_with(partner)
                
    def action(self, lattice):
        '''
        Runs all actions a monomer can perform in succession.
        '''
        self.diffuse(lattice)
        self.rotate(lattice)
        self.couple(lattice)