Simulating crystalline LiH with OpenFermion
===========================================

This project demonstrates how `OpenFermion <https://github.com/quantumlib/OpenFermion>`__ and `ProjectQ <https://projectq.ch/>`__ can be used to calculate the energy of crystalline solids using the variational quantum eigensolver (VQE) algorithm for quantum computers.

The fermionic Hamiltonian is constructed in the dual basis and mapped to a qubit Hamiltonian using the Jordan-Wigner transformation.
To reduce the computational complexity, the active space of the Hamiltonian is reduced by freezing 14 of the 16 electrons present in the 8-atom conventional unit cell of LiH.
The energy is then optimized using the VQE with a UCCSD ansatz, using the ProjectQ backend to simulate the quantum circuits.

The graph below shows how the ground-state energy decreases as the number of orbitals included in the coupled-cluster ansatz is increased.
(Note that two ative orbitals corresponds to a wavefunction comprised of a single Slater determinant.)
The results show that the ground-state energy is not well converged with eight orbitals in the basis. 
Although the calculations are very coarse (in terms of basis set, active space, and supercell size), they demonstrate the basic workflow for calculating energies of crystalline solids.

.. image:: https://github.com/max-radin/LiH-example/raw/master/energies.png
    :align: left
    
