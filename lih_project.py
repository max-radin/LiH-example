from openfermion.hamiltonians import plane_wave_hamiltonian
from openfermion.utils import Grid
from openfermion.transforms import jordan_wigner
import os
from scipy.optimize import minimize

from openfermion.config import *
from openfermionprojectq import *

from openfermion.utils import uccsd_singlet_paramsize
from openfermion.utils import count_qubits
from projectq.ops import X, All, Measure
from projectq.backends import CommandPrinter, CircuitDrawer

class LiH_Project:
    """Class describing a VQE/UCCSD calculation on crystalline LiH.
    
    Attributes:
        a (float): lattice parameter (Bohr)
        n (int): number of grid subdivisions in each direction
        n_active_el (int): number of active electrons
        n_active_orb (int): number of active orbitals
        N_units (int): number of formula units in the cell
        opt_energy (float): optimized energy
        opt_amplitudes (list): optimized UCCSD amplitudes
    """
    def __init__(self, a=7.72, n=3, n_active_el=2, n_active_orb=4):
        self.a = a
        self.n = n
        self.n_active_el = n_active_el
        self.n_active_orb = n_active_orb
        self.N_units = 4
        self.opt_energy = None
        self.opt_amplitudes = None

        species_a = 'Li'
        species_b = 'H'

        # Construct a fermion Hamiltonian
        grid = Grid(3, self.n, self.a)
        geometry = [(species_a, (0, 0, 0)),
                    (species_a, (0, 0.5, 0.5)),
                    (species_a, (0.5, 0, 0.5)),
                    (species_a, (0.5, 0.5, 0)),
                    (species_b, (0.5, 0.5, 0.5)),
                    (species_b, (0.5, 0, 0)),
                    (species_b, (0, 0.5, 0)),
                    (species_b, (0, 0, 0.5))
                    ]
        self.fermion_hamiltonian = plane_wave_hamiltonian(
                grid,
                geometry=geometry,
                spinless=False,
                plane_wave=False,
                include_constant=False,
                e_cutoff=None)

        # Freeze specified orbitals
        n_qubits = 2*n*n*n
        n_electrons = 4*3+4*1
        self.fermion_hamiltonian.freeze_orbitals(
                range(n_electrons-self.n_active_el),
                range(n_electrons-self.n_active_el+self.n_active_orb,
                      n_qubits))

        # Construct qubit Hamiltonian
        self.qubit_hamiltonian = jordan_wigner(self.fermion_hamiltonian)

        # Initialize com;iler engine
        self.compiler_engine = uccsd_trotter_engine()

    def energy_objective(self, packed_amplitudes):
        """Evaluate the energy of a UCCSD singlet wavefunction with packed_amplitudes
        Args:
            packed_amplitudes(ndarray): Compact array that stores the unique
                amplitudes for a UCCSD singlet wavefunction.

        Returns:
            energy(float): Energy corresponding to the given amplitudes
        """
        os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

        # Set Jordan-Wigner initial state with correct number of electrons
        wavefunction = self.compiler_engine.allocate_qureg(self.n_active_orb)
        for i in range(self.n_active_el):
            X | wavefunction[i]

        # Build the circuit and act it on the wavefunction
        evolution_operator = uccsd_singlet_evolution(packed_amplitudes, 
                                                     self.n_active_orb,
                                                     self.n_active_el)
        evolution_operator | wavefunction
        self.compiler_engine.flush()

        # Evaluate the energy and reset wavefunction
        energy = self.compiler_engine.backend.get_expectation_value(
                self.qubit_hamiltonian, 
                wavefunction)
        All(Measure) | wavefunction
        self.compiler_engine.flush()
        return energy
    
    def get_energy(self):
        """Minimize the energy objective function."""
        
        n_amplitudes = uccsd_singlet_paramsize(self.n_active_orb,
                                               self.n_active_el)
        initial_amplitudes = [0.001 for i in range(n_amplitudes)]

        # Run VQE Optimization to find new CCSD parameters
        opt_result = minimize(self.energy_objective, initial_amplitudes,
                              method="CG", options={'disp':True})
        self.opt_energy, self.opt_amplitudes = opt_result.fun, opt_result.x
        return self.opt_energy/self.N_units
