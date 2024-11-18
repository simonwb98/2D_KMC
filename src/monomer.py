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
        self.coupling_energy = float(coupling_energy)
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
        x, y = self.get_position()
        neighbors = lattice.get_cached_neighbors(x, y, self.get_orientation())
        unoccupied_sites = [site for site in neighbors if not lattice.is_occupied(*site)]
        if unoccupied_sites:
            x_new, y_new = random.choice(unoccupied_sites)
            lattice.move_monomer(self, x_new, y_new)


    def rotate(self, lattice):
        rotation_prob = self.rotation_probability(lattice)

        if random.random() < rotation_prob:
            self.set_orientation(random.choice([o for o in self.orientations if not o == self.orientation]))

    def couple(self, lattice):
        """
        Perform coupling using cached next-nearest neighbors.
        """
        x, y = self.get_position()
        next_nearest = lattice.get_cached_next_nearest_neighbors(x, y, self.orientation)
        candidates = [
            lattice.grid[ny][nx] for nx, ny in next_nearest
            if lattice.grid[ny][nx] and lattice.grid[ny][nx].orientation != self.orientation
        ]
        if candidates:
            partner = random.choice(candidates)
            self.couple_with(partner)

            
    def calculate_diffusion_rate(self, lattice):
        """
        Calculate the total diffusion rate for this monomer based on unoccupied neighbors.

        Args:
            lattice (Lattice): The lattice object.

        Returns:
            float: Total diffusion rate.
        """
        if not self.coupled:
            temperature = lattice.temperature
            base_rate = self.diffusion_rate * math.exp(-self.diffusion_energy / (k_B * temperature))
            neighbors = lattice.get_neighbours(*self.get_position())
            unoccupied_sites = [site for site in neighbors if not lattice.is_occupied(*site)]
            return len(unoccupied_sites) * base_rate
        else: 
            return 0

    def calculate_coupling_rate(self, lattice):
        """
        Calculate the total coupling rate for this monomer based on valid next-nearest neighbors.

        Args:
            lattice (Lattice): The lattice object.

        Returns:
            float: Total coupling rate.
        """
        temperature = lattice.temperature
        base_rate = self.coupling_rate * math.exp(-self.coupling_energy / (k_B * temperature))

        neighbours = lattice.get_neighbours(*self.get_position())
        if any(lattice.is_occupied(nx, ny) for (nx, ny) in neighbours):
            return 0 # disallow coupling when there are nearest neighbours to this monomer. This condition basically realizes the fact that monomers physically restrict each other (geometric hindrance) - a cleaner way of doing this is to disallow diffusion into sites that have monomers that are nearest neighbours, but this is fine also.
        next_neighbours = lattice.get_next_nearest_neighbours(*self.get_position(), self.get_orientation()) # this could potentially be made faster sometime down the line
        valid_partners = [lattice.grid[ny][nx] for (nx, ny) in next_neighbours if not lattice.grid[ny][nx] == None and lattice.grid[ny][nx].get_orientation() != self.get_orientation()]
        return len(valid_partners) * base_rate

    def calculate_rotation_rate(self, lattice):
        """
        Calculate the rotation rate for this monomer.
        Currently independent of the local environment.

        Args:
            lattice (Lattice): The lattice object.

        Returns:
            float: Rotation rate.
        """
        temperature = lattice.temperature
        return 2 * self.rotation_rate * math.exp(-self.rotation_energy / (k_B * temperature)) if not self.coupled else 0

    def calculate_total_rate(self, lattice):
        """
        Calculate the total rate for this monomer by summing rates for all actions.

        Args:
            lattice (Lattice): The lattice object.

        Returns:
            float: Total rate for the monomer.
        """
        return (
            self.calculate_diffusion_rate(lattice) +
            self.calculate_coupling_rate(lattice) +
            self.calculate_rotation_rate(lattice)
        )


    def update_rates(self, lattice):
        """
        Update the rates for diffusion, rotation, and coupling based on the local environment.
        Skip updating if the monomer is already coupled.
        """
        if self.coupled:
            self.cached_diffusion_rate = 0
            self.cached_rotation_rate = 0
            self.cached_coupling_rate = 0
        else:
            self.cached_diffusion_rate = self.calculate_diffusion_rate(lattice)
            self.cached_rotation_rate = self.calculate_rotation_rate(lattice)
            self.cached_coupling_rate = self.calculate_coupling_rate(lattice)


    def calculate_total_rate(self):
        """
        Use cached rates to compute the total rate for the monomer.
        """
        return self.cached_diffusion_rate + self.cached_rotation_rate + self.cached_coupling_rate
