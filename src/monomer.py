# src/monomer.py

import math, random
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

class Monomer:
    def __init__(self, monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, orientation = 0):
        self.monomer_type = monomer_type
        self.diffusion_rate = diffusion_rate          
        self.diffusion_energy = diffusion_energy
        self.rotation_rate = rotation_rate
        self.rotation_energy = rotation_energy
        self.coupled = False
        self.position = None

    def set_position(self, x, y):
        self.position = (x, y)

    def get_position(self):
        return self.position
    
    def diffusion_probability(self, lattice):
        temperature = lattice.temperature
        return self.diffusion_rate * math.exp(-self.diffusion_energy / (k_B * temperature))
    
    def rotation_probability(self, lattice):
        temperature = lattice.temperature
        return self.rotation_rate * math.exp(-self.rotation_energy / (k_B * temperature))
    
    def diffuse(self, lattice):
        diffusion_prob = self.diffusion_probability(lattice)
        neighbours = lattice.get_neighbours(lattice, *self.get_position)

        if random.random() < diffusion_prob:
            x_new, y_new = random.choice(neighbours)
            self.set_position(x_new, y_new)
