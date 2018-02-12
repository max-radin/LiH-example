#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 15:22:52 2018

@author: max
"""

import openfermion as of

a = 4.085
species_a = 'Li'
species_b = 'H'

grid = of.utils.Grid(3, 3, a)

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
                                           include_constant=True,
                                           e_cutoff=None)

print('Hamiltonian has {} terms'.format(len(H.terms)))

H_jw = of.transforms.jordan_wigner(H)

print('Transformated Hamiltonian has {} terms'.format(len(H_jw.terms)))
