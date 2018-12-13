# -*- coding: utf-8 -*-
"""
.. module:: geek
   :platform: Unix, Windows
   :synopsis: GEneralised Elementary Kinetics

.. moduleauthor:: geek team

[---------]

Copyright 2018 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""
from ecell4 import *
from ecell4.bd import BDWorld, BDSimulator
from ecell4.egfrd import EGFRDWorld, EGFRDSimulator
from ecell4.core import GSLRandomNumberGenerator
import numpy as np
import random

AVOGADRO_NUMBER = 6e23

"""
The following script shall provide a validation of GEEK and the brownian reaction dynamics by comparing its results with 
other approaches for uni-molecular and for bi-molecular reaction

The simulations are conduced using a single sized crowding size for an efficient comparison with the ecell4 BRD
and the EGFRD 
"""


"""
Diffusion limited conditions
"""

# gamma = 4*pi*4e-9*200e-12*1000*6e23 ~ 6e9

parameters_diff_lim = {
    'K_eq':  50e-6,     # M
    'k_fwd': 5e9,       # 1/Ms
    'r_A': 2e-9,        # m
    'r_B': 2e-9,        # m
    'r_C': 3e-9,        # m
    'D_A': 100e-12,     # m^2/s
    'D_B': 100e-12,     # m^2/s
    'D_C':  75e-12,     # m^2/s
    'm_A': 10,          # kDa
    'm_B': 10,          # kDa
    'm_C': 20,          # kDa
    'A_0': 50e-6,       # M
    'B_0': 50e-6,       # M
    'C_0': 0,           # M
    'volume': 10e-18,   # L
    't_max': 1e-5,      # s
    'dt': 0.25e-9,      # s
}




"""
Reaction limited conditions
"""

# gamma = 4*pi*4e-9*200e-12*1000*6e23 ~ 6e9

parameters_reaction_lim = {
    'K_eq':  50e-6,     # M
    'k_fwd': 1e7,       # 1/Ms
    'r_A': 2e-9,        # m
    'r_B': 2e-9,        # m
    'r_C': 3e-9,        # m
    'D_A': 100e-12,     # m^2/s
    'D_B': 100e-12,     # m^2/s
    'D_C':  75e-12,     # m^2/s
    'm_A': 10,          # kDa
    'm_B': 10,          # kDa
    'm_C': 20,          # kDa
    'A_0': 50e-6,       # M
    'B_0': 50e-6,       # M
    'C_0': 0,           # M
    'volume': 10e-18,   # L
    't_max': 1e-6,      # s
    'dt': 0.25e-9,      # s
    'mu_mass': 31.9     # s
}


""" 
Helper functions 
"""
def mass2rad(mass):
   radius = 0.0515*(mass*1000)**(0.393) # Mass in kDa
   return radius


def rad2mass(radius):
   M = (radius/0.0515)**(1./0.393)/1000.0 #Radius in nm
   return M

def rad2diff(radius):
    viscosity = 0.7e-3 # Pa s
    temperatur = 310.15 # K
    kbT = temperatur*1.38064852e-23
    D = kbT/(6*np.pi*viscosity*radius) #Radius in m
    return D # in m^2/s




""" 
Simulation functions
"""

def geek_simulations(parameters,phi= 0.0):
    pass


def openbread_simulation(parameters,phi= 0.0):
    pass


def ecell4_gfrd_simulation(parameters,phi= 0.0, seed=1):
    m = NetworkModel()

    rng = GSLRandomNumberGenerator()
    rng.seed(seed)

    #Rescale to mum
    s = 1e6
    s2 = 1e12
    s3 = 1e18

    A = Species('A', parameters['r_A']*s, parameters['D_A']*s2)
    B = Species('B', parameters['r_B']*s, parameters['D_B']*s2)
    C = Species('C', parameters['r_C']*s, parameters['D_C']*s2)

    if phi > 0:
        R_crw = mass2rad(parameters['mu_mass'])
        volume_crw = 4.0 / 3.0 * np.pi * R_crw ** 3
        crw = Species('crw',R_crw*s, rad2diff(R_crw)*s2 )

    gamma = 4.0*np.pi*(parameters['r_A']+parameters['r_B'])*(parameters['D_A']+parameters['D_B'])*s3
    rescaled_keff = parameters['k_fwd']/AVOGADRO_NUMBER/1000.0*s3

    effective_k_binding = gamma*rescaled_keff/(gamma-rescaled_keff)

    print(gamma)
    print(effective_k_binding)
    m.add_reaction_rule(create_binding_reaction_rule(A, B, C, effective_k_binding))

    m.add_reaction_rule(create_unbinding_reaction_rule(C, A, B, parameters['k_fwd']*parameters['K_eq']))

    a = (parameters['volume']/1000)**(1/3) * s # in mum

    w=EGFRDWorld(edge_lengths=Real3(a,a,a)) #,rng=rng)

    w.bind_to(m)

    # Add the species in the concentrations
    N_A = parameters['A_0']*parameters['volume']*AVOGADRO_NUMBER
    N_B = parameters['B_0']*parameters['volume']*AVOGADRO_NUMBER
    N_C = parameters['C_0']*parameters['volume']*AVOGADRO_NUMBER

    w.add_molecules(A, N_A)
    w.add_molecules(B, N_B)
    w.add_molecules(C, N_C)


    # Add crowding species
    if phi > 0:
        N_crw = round(parameters['volume']*phi/volume_crw)
        w.add_molecules(crw, N_crw)

    obs = FixedIntervalNumberObserver(parameters['dt']*10.0, ['A', 'B', 'C'])
    sim = EGFRDSimulator(w)
    sim.run(parameters['t_max'], obs)
    return obs.data()


def ecell4_brd_simulation(parameters,phi=0.0,seed=1):
    m = NetworkModel()

    rng = GSLRandomNumberGenerator()
    rng.seed(seed)

    #Rescale to mum
    s = 1e6
    s2 = 1e12
    s3 = 1e18
    A = Species('A', parameters['r_A']*s, parameters['D_A']*s2)
    B = Species('B', parameters['r_B']*s, parameters['D_B']*s2)
    C = Species('C', parameters['r_C']*s, parameters['D_C']*s2)

    if phi > 0:
        R_crw = mass2rad(parameters['mu_mass'])
        volume_crw = 4.0 / 3.0 * np.pi * R_crw ** 3
        crw = Species('crw',R_crw*s, rad2diff(R_crw)*s2 )

    gamma = 4.0*np.pi*(parameters['r_A']+parameters['r_B'])*(parameters['D_A']+parameters['D_B'])*s3
    rescaled_keff = parameters['k_fwd']/AVOGADRO_NUMBER/1000.0*s3
    effective_k_binding = gamma*rescaled_keff/(gamma-rescaled_keff)

    print(gamma)
    print(effective_k_binding)
    m.add_reaction_rule(create_binding_reaction_rule(A, B, C, effective_k_binding))

    m.add_reaction_rule(create_unbinding_reaction_rule(C, A, B, parameters['k_fwd']*parameters['K_eq']))

    m = get_model()

    a = (parameters['volume'] / 1000) ** (1 / 3) * s  # in m
    w = BDWorld(Real3(a, a, a), rng=rng)
    w.bind_to(m)

    # Add the species in the concentrations
    N_A = parameters['A_0']*parameters['volume']*AVOGADRO_NUMBER
    N_B = parameters['B_0']*parameters['volume']*AVOGADRO_NUMBER
    N_C = parameters['C_0']*parameters['volume']*AVOGADRO_NUMBER


    w.add_molecules(A, N_A)
    w.add_molecules(B, N_B)
    w.add_molecules(C, N_C)

    # Add crowding species
    if phi > 0:
        N_crw = round(parameters['volume']*phi/volume_crw)
        w.add_molecules(crw, N_crw)

    obs = FixedIntervalNumberObserver(parameters['dt']*10.0, ['A', 'B', 'C'])
    sim = BDSimulator(w)
    #sim.set_dt(parameters['dt'])
    sim.run(parameters['t_max'], obs)
    return obs.data()





