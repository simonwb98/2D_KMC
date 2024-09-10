# src/polymer.py

class Polymer:
    def __init__(self, monomers):
        self.monomers = monomers
        self.occupied_cells = [monomer.position for monomer in monomers]

    def add_monomer(self, monomer):
        if monomer not in self.monomers:
            self.monomers.append(monomer)
            self.occupied_cells.append(monomer.position)
            monomer.coupled = True

    def grow(self, lattice):
        unique_neighbours = set()

        for monomer in self.monomers:
            neighbours = lattice.get_neighbours(*monomer.position)
            unique_neighbours.update(neighbours)

        for (nx, ny) in unique_neighbours:
            if lattice.is_occupied(nx, ny) and not (nx, ny) in self.occupied_cells: # check if the unique neighbours are (a) occupied, and (b) not part of the existing network
                self.add_monomer(lattice.grid[ny][nx]) # add monomer instance to the polymer
                self.occupied_cells.append(lattice.grid[ny][nx].position)

            