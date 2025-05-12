# src/monomer.py

import math, random
k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

class Monomer:
    def __init__(self, monomer_type, diffusion_rate, diffusion_energy, rotation_rate, rotation_energy, coupling_rate, coupling_energy, dehalogen_rate, dehalogen_energy, orientations = [0, 180]):
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
        self.site1 = True
        self.site2 = True
        self.site3 = True
        self.rotations = 0
        self.dehalogen_energy = dehalogen_energy
        self.dehalogen_rate = dehalogen_rate

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
        neighbors = lattice.get_neighbours(x, y)
        unoccupied_sites = [site for site in neighbors if not lattice.is_occupied(*site)]
        if unoccupied_sites:
            x_new, y_new = random.choice(unoccupied_sites)
            lattice.move_monomer(self, x_new, y_new)

    def rotate(self, lattice):
        rotation_prob = self.calculate_rotation_rate(lattice)
#        print("id like to ROTATE")
        if random.random() < rotation_prob:
            self.set_orientation(random.choice([o for o in self.orientations if not o == self.orientation]))
            self.rotations += 1
#            print("ROTATE")
    def get_halogenation(self, lattice, partner):
        
        """This function checks the halogenation of the monomers trying to couple. Returns 0 if coupling is forbidden, 1 if it is allowed."""
        partner_neighbours = lattice.get_neighbours(*partner.get_position())
        orientation1 = self.rotations % 6
        x1, y1 = self.get_position()
        if y1%2:
            x1 = x1 + 0.5 # shifting the x coordinate for the angle calculation
        site1_1 = self.site1        
        site2_1 = self.site2
        site3_1 = self.site3
        
        orientation2 = partner.rotations % 6
        x2, y2 = partner.get_position()
        if y2%2:
            x2 = x2 + 0.5 # shifting the x coordinate for the angle calculation
        site1_2 = partner.site1
        site2_2 = partner.site2
        site3_2 = partner.site3
        
        angle  = int(math.degrees(math.atan((x2-x1)/(y2-y1))))
#        print(orientation1, orientation2, angle)        
        
        # now use the orientation and the angle to determine the sites involved in the bonding
        if (orientation1 == 0):
            if (angle == 0):
                if (orientation2 == 1):
                    return site1_1*site2_2                
                if (orientation2 == 3):
                    return site1_1*site1_2
                if (orientation2 == 5):
                    return site1_1*site3_2
            if (angle == 56):
                if (orientation2 == 1):
                    return site2_1*site3_2    
                if (orientation2 == 3):
                    return site2_1*site2_2
                if (orientation2 == 5):
                    return site2_1*site1_2
            if (angle == -56):
                if (orientation2 == 1):
                    return site3_1*site1_2    
                if (orientation2 == 3):
                    return site3_1*site3_2
                if (orientation2 == 5):
                    return site3_1*site2_2

        if (orientation1 == 1):
            if (angle == 0):
                if (orientation2 == 0):
                    return site2_1*site1_2
                if (orientation2 == 2):
                    return site2_1*site3_2
                if (orientation2 == 4):
                    return site2_1*site2_2
            if (angle == 56):
                if (orientation2 == 0):
                    return site1_1*site3_2    
                if (orientation2 == 2):
                    return site1_1*site2_2
                if (orientation2 == 4):
                    return site1_1*site1_2
            if (angle == -56):
                if (orientation2 == 0):
                    return site3_1*site2_2
                if (orientation2 == 2):
                    return site3_1*site1_2
                if (orientation2 == 4):
                    return site3_1*site3_2

        if (orientation1 == 2):
            if (angle == 0):
                if (orientation2 == 1):
                    return site3_1*site2_2
                if (orientation2 == 3):
                    return site3_1*site1_2
                if (orientation2 == 5):
                    return site3_1*site3_2
            if (angle == 56):
                if (orientation2 == 1):
                    return site1_1*site3_2
                if (orientation2 == 3):
                    return site1_1*site2_2
                if (orientation2 == 5):
                    return site1_1*site1_2
            if (angle == -56):
                if (orientation2 == 1):
                    return site2_1*site1_2
                if (orientation2 == 3):
                    return site2_1*site3_2
                if (orientation2 == 5):
                    return site2_1*site2_2

        if (orientation1 == 3):
            if (angle == 0):
                if (orientation2 == 0):
                    return site1_1*site1_2
                if (orientation2 == 2):
                    return site1_1*site3_2
                if (orientation2 == 4):
                    return site1_1*site2_2
            if (angle == 56):
                if (orientation2 == 0):
                    return site3_1*site3_2
                if (orientation2 == 2):
                    return site3_1*site2_2
                if (orientation2 == 4):
                    return site3_1*site1_2
            if (angle == -56):
                if (orientation2 == 0):
                    return site2_1*site2_2
                if (orientation2 == 2):
                    return site2_1*site1_2
                if (orientation2 == 4):
                    return site2_1*site3_2

        if (orientation1 == 4):
            if (angle == 0):
                if (orientation2 == 1):
                    return site2_1*site2_2
                if (orientation2 == 3):
                    return site2_1*site1_2
                if (orientation2 == 5):
                    return site2_1*site3_2
            if (angle == 56):
                if (orientation2 == 1):
                    return site3_1*site3_2
                if (orientation2 == 3):
                    return site3_1*site2_2
                if (orientation2 == 5):
                    return site3_1*site1_2
            if (angle == -56):
                if (orientation2 == 1):
                    return site1_1*site1_2
                if (orientation2 == 3):
                    return site1_1*site3_2
                if (orientation2 == 5):
                    return site1_1*site2_2

        if (orientation1 == 5):
            if (angle == 0):
                if (orientation2 == 0):
                    return site1_1*site1_2
                if (orientation2 == 2):
                    return site1_1*site3_2
                if (orientation2 == 4):
                    return site1_1*site2_2
            if (angle == 56):     
                if (orientation2 == 0):
                    return site2_1*site3_2
                if (orientation2 == 2):
                    return site2_1*site2_2
                if (orientation2 == 4):
                    return site2_1*site1_2
            if (angle == -56):
                if (orientation2 == 0):
                    return site3_1*site2_2
                if (orientation2 == 2):
                    return site3_1*site1_2
                if (orientation2 == 4):
                    return site3_1*site3_2
        
    def calculate_dehalogen_rate(self, lattice):
        temperature = lattice.temperature
        rate = self.dehalogen_rate*math.exp(-self.dehalogen_energy/(k_B * temperature))
        return rate

    def dehalogenate(self, lattice):
        '''Checks if each halogenation site should be dehalogenated'''
        rate = self.calculate_dehalogen_rate(lattice)
        sites = [self.site1, self.site2, self.site3]
        for i in range(len(sites)):
            if not sites[i]:
                continue

            if random.random() < rate:
                sites[i] = False
        self.site1=sites[0]
        self.site2=sites[1]
        self.site3=sites[2]

        return

    def couple(self, lattice):
        """
        Perform coupling using cached next-nearest neighbors.
        """
        c_rate = self.calculate_coupling_rate(lattice)
        if random.random() < c_rate:
#            print("INITIATE COUPLE")
            x, y = self.get_position()
            next_nearest = lattice.get_next_nearest_neighbours(x, y, self.orientation)
            candidates = [
                lattice.grid[ny][nx] for nx, ny in next_nearest
                if lattice.grid[ny][nx] and lattice.grid[ny][nx].orientation != self.orientation
            ]
            if candidates:
                partner = random.choice(candidates)
                halogen_bool = self.get_halogenation(lattice, partner)
#                print(halogen_bool)
                if halogen_bool==1 or halogen_bool==None:
                    return
                self.couple_with(partner)
                print("COUPLED")

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
        if self.coupled: 
            return 0
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
        if self.coupled:
            return 0 
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
