import openfermion as of


def build_hamiltonian(a, n):
    # Construct a Fermion Hamiltonian for crystalline LiH.
    #
    # Args:
    #   a: float indicating the lattice parameter
    #   n: number of grid subdivisions in each direction

    species_a = 'Li'
    species_b = 'H'

    grid = of.utils.Grid(3, n, a)

    geometry = [(species_a, (0, 0, 0)),
                (species_a, (0, 0.5, 0.5)),
                (species_a, (0.5, 0, 0.5)),
                (species_a, (0.5, 0.5, 0)),
                (species_b, (0.5, 0.5, 0.5)),
                (species_b, (0.5, 0, 0)),
                (species_b, (0, 0.5, 0)),
                (species_b, (0, 0, 0.5))
                ]

    H = of.hamiltonians.plane_wave_hamiltonian(grid,
                                               geometry=geometry,
                                               spinless=False,
                                               plane_wave=False,
                                               include_constant=False,
                                               e_cutoff=None)

    return(H)
