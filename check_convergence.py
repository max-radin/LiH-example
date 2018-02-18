from lih_project import LiH_Project
import matplotlib.pyplot as plt

# Calculate the energy of crystalline LiH, freezing all but two electrons.

# %% Calculate energy

n_active_orb_range = [2, 4, 6, 8]

projs = []  # Calculation projects
energies = []  # Optimized energies

# Calculate the energy for a single Fock state explcitly
projs.append(LiH_Project(n_active_el=0, n_active_orb=0))
energies.append(projs[0].fermion_hamiltonian.terms[()].real/4.0)

# Use the VQE to calculate the energy more accurately
for n_active_orb in n_active_orb_range[1:]:
    projs.append(LiH_Project(n_active_orb=n_active_orb))
    energies.append(projs[-1].get_energy())

# %% Plot results
plt.figure()
plt.plot(n_active_orb_range, energies, 'o-')
plt.xlabel('# of active orbitals')
plt.ylabel('Energy (Hartree/LiH)')
plt.show()
plt.savefig('energies.png')
