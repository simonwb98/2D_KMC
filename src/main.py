from lattice import Lattice

def main():
    lattice = Lattice(width = 3, rotational_symmetry = 6, periodic=True)

    x, y = 2, 2
    neighbours = lattice.get_neighbours(x, y)
    print(f"Neighbours of ({x, y}): {neighbours}")


if __name__ == "__main__":
    main()