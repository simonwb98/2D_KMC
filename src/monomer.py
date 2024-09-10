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
        self.energy = coupling_energy
        self.coupled = False
        self.position = None
        self.orientations = orientations
        self.orientation = random.choice(orientations)

    def set_position(self, x, y):
        self.position = (x, y)

    def get_position(self):
        return self.position
    
    def diffusion_probability(self, lattice):
        temperature = lattice.temperature
        return self.diffusion_rate * math.exp(-self.diffusion_energy / (k_B * temperature)) if not self.coupled else 0 # for two or more monomers, the diffusion probability is (for now) set to 0.
    
    def rotation_probability(self, lattice):
        temperature = lattice.temperature
        return self.rotation_rate * math.exp(-self.rotation_energy / (k_B * temperature)) if not self.coupled else 0
    
    def coupling_probability(self, lattice):
        temperature = lattice.temperature
        return self.coupling_rate * math.exp(-self.coupling_energy / k_B * temperature) if not self.coupled else 0
    
    def diffuse(self, lattice):
        diffusion_prob = self.diffusion_probability(lattice)
        neighbours = lattice.get_neighbours(*self.get_position())

        if random.random() < diffusion_prob:
            x_new, y_new = random.choice(neighbours)
            lattice.move_monomer(self, x_new, y_new)

    def rotate(self, lattice):
        rotation_prob = self.rotation_probability(lattice)

        if random.random() < rotation_prob:
            new_orientation = random.choice([o for o in self.orientations if not o == self.orientation])


    def couple(self, lattice):
        coupling_prob = self.coupling_probability(lattice)
