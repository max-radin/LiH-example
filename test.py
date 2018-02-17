import os

from scipy.optimize import minimize

from openfermion.config import *
from openfermionprojectq import *

from openfermion.transforms import jordan_wigner
from openfermion.utils import uccsd_singlet_paramsize
from openfermion.utils import count_qubits
from projectq.ops import X, All, Measure
from projectq.backends import CommandPrinter, CircuitDrawer

from LiH_utils import build_hamiltonian

# %% Construct qubit Hamiltonian
n = 3  # Number of grid subdivisions in each dimension
a = 7.72 # Lattice parameter (Bohr)
n_qubits = 2*n*n*n
n_electrons = 4*3+4*1

fermion_hamiltonian = build_hamiltonian(a, n)

print('Fermion Hamiltonian has {} terms'.format(
        len(fermion_hamiltonian.terms)))

n_active_el = 2  # Number of active electrons
n_active_orb = 4  # Number of active orbitals
fermion_hamiltonian.freeze_orbitals(range(n_electrons-n_active_el),
                                    range(n_electrons-n_active_el+n_active_orb,
                                          n_qubits))

print('Frozen Fermion Hamiltonian has {} terms'.format(
        len(fermion_hamiltonian.terms)))

ref_energy = fermion_hamiltonian.terms[()]
print("Reference energy: {}".format(ref_energy))

qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)

compiler_engine = uccsd_trotter_engine()


print('Qubit Hamiltonian has {} terms'.format(
        len(qubit_hamiltonian.terms)))
print("Qubit Hamiltonian acts on {} qubits".format(
        count_qubits(qubit_hamiltonian)))


# %% Define objective function
def energy_objective(packed_amplitudes):
    """Evaluate the energy of a UCCSD singlet wavefunction with packed_amplitudes
    Args:
        packed_amplitudes(ndarray): Compact array that stores the unique
            amplitudes for a UCCSD singlet wavefunction.

    Returns:
        energy(float): Energy corresponding to the given amplitudes
    """
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Set Jordan-Wigner initial state with correct number of electrons
    wavefunction = compiler_engine.allocate_qureg(n_active_orb)
    for i in range(n_active_el):
        X | wavefunction[i]

    # Build the circuit and act it on the wavefunction
    evolution_operator = uccsd_singlet_evolution(packed_amplitudes, 
                                                 n_active_orb,
                                                 n_active_el)
    evolution_operator | wavefunction
    compiler_engine.flush()

    # Evaluate the energy and reset wavefunction
    energy = compiler_engine.backend.get_expectation_value(qubit_hamiltonian, 
                                                           wavefunction)
    All(Measure) | wavefunction
    compiler_engine.flush()
    return energy

# %% Minimize the objective function
n_amplitudes = uccsd_singlet_paramsize(n_active_orb, n_active_el)
initial_amplitudes = [0.001 for i in range(n_amplitudes)]
initial_energy = energy_objective(initial_amplitudes)

# Run VQE Optimization to find new CCSD parameters
opt_result = minimize(energy_objective, initial_amplitudes,
                      method="CG", options={'disp':True})

opt_energy, opt_amplitudes = opt_result.fun, opt_result.x
print("\nOptimal UCCSD Singlet Energy: {}".format(opt_energy+ref_energy))
print("Optimal UCCSD Singlet Amplitudes: {}".format(opt_amplitudes))
print("Initial Energy of UCCSD with CCSD amplitudes: {} Hartrees".
      format(initial_energy+ref_energy))
