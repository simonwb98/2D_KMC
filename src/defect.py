import math, random
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

class Defect:
    def __init__(self, diffusion_rate, diffusion_energy, nucleation_prob):
        self.diffusion_rate = diffusion_rate          
        self.diffusion_energy = diffusion_energy
        self.nucleation_prob = nucleation_prob
        self.nucleating = False
        self.position = None
        self.coupled = False

    def set_position(self, x, y):
        self.position = (x, y)

    def get_orientation(self):
        return False
    
    def get_position(self):
        return self.position
    
    def diffusion_probability(self, lattice):
        temperature = lattice.temperature
        return self.diffusion_rate * math.exp(-self.diffusion_energy / (k_B * temperature)) if not self.nucleating else 0 # for two or more coupled monomers, the diffusion probability is (for now) set to 0.
    
    def couple_with(self, other):
        self.coupled = True
        other.nucleating = True

    def diffuse(self, lattice):
        diffusion_prob = self.diffusion_probability(lattice) # get probability
        
        if random.random() < diffusion_prob and not self.nucleating: # based on the probability, decide if diffuse or not
            neighbours = lattice.get_neighbours(*self.get_position())
            x_new, y_new = random.choice(neighbours)
            lattice.move_defect(self, x_new, y_new)
    
    
    def action(self, lattice):
        if not self.nucleating:
            self.diffuse(lattice)
           
        else: pass