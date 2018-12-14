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


import numpy as np
import random

from pandas import DataFrame

AVOGADRO_NUMBER = 6e23

"""
The following script shall provide a validation of GEEK and the brownian reaction dynamics by comparing its results with 
other approaches for uni-molecular and for bi-molecular reaction

The simulations are conduced using a single sized crowding size for an efficient comparison with the ecell4 BRD
and the EGFRD 
"""


"""
Diffusion controlled conditions
"""

# gamma = 4*pi*4e-9*200e-12*1000*6e23 ~ 6e9

parameters_diff_lim = {
    'K_eq':  50e-6,     # M
    'k_fwd': 5e9,       # 1/Ms
    'r_A': 2e-9,        # m
    'r_B': 2e-9,        # m
    'r_C': 3e-9,        # m
    'D_A': 500e-12,     # m^2/s
    'D_B': 500e-12,     # m^2/s
    'D_C': 350e-12,     # m^2/s
    'm_A': 10,          # kDa
    'm_B': 10,          # kDa
    'm_C': 20,          # kDa
    'A_0': 50e-6,       # M
    'B_0': 50e-6,       # M
    'C_0': 0,           # M
    'volume': 10e-18,   # L
    't_max': 1e-5,      # s
    'dt': 0.25e-9,      # s
    'mu_mass': 21.1     # s
}




"""
Reaction limited conditions
"""

# gamma = 4*pi*4e-9*200e-12*1000*6e23 ~ 6e9

parameters_reaction_lim = {
    'K_eq':  50e-6,     # M
    'k_fwd': 5e7,       # 1/Ms
    'r_A': 2e-9,        # m
    'r_B': 2e-9,        # m
    'r_C': 3e-9,        # m
    'D_A': 500e-12,     # m^2/s
    'D_B': 500e-12,     # m^2/s
    'D_C': 350e-12,     # m^2/s
    'm_A': 10,          # kDa
    'm_B': 10,          # kDa
    'm_C': 20,          # kDa
    'A_0': 50e-6,       # M
    'B_0': 50e-6,       # M
    'C_0': 0,           # M
    'volume': 10e-18,   # L
    't_max': 1e-3,      # s
    'dt': 0.25e-9,      # s
    'mu_mass': 21.1     # s
}


""" 
Helper functions 
"""
def mass2rad(mass):
   radius = 0.0515*(mass*1000)**(0.393) # Mass in kDa
   return radius #Radius in nm


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


def openbread_simulation(parameters, phi= 0.0, seed=1):

    from openbread.core import Species,ParticleModel,Reaction
    # Construct species in the model

    A = Species(name='A',
                diffusion_constant=parameters['D_A'],
                collision_radius=parameters['r_A'],
                mass=parameters['m_A'],)
    B = Species(name='B',
                diffusion_constant=parameters['D_B'],
                collision_radius=parameters['r_B'],
                mass=parameters['m_B'],)

    C = Species(name='C',
                diffusion_constant=parameters['D_C'],
                collision_radius=parameters['r_C'],
                mass=parameters['m_C'],)

    species = [A, B, C, ]

    # Define microscopic reaction rate constants:
    k1f = parameters['k_fwd']  # 1/Ms
    k1b = parameters['k_fwd'] * parameters['K_eq']  # 1/s

    # Setup particle simulation environemnt
    volume = parameters['volume'] # (0.1 mum)^3 in L

    medium = ParticleModel.Medium(viscosity=0.7e-3,  # Pa s
                                  temperatur=310.15)

    crowding = ParticleModel.Crowding(volume_fraction=phi,
                                      mu=np.log(parameters['mu_mass']),
                                      sigma=0,
                                      max_size=10e-3)

    particle_model = ParticleModel(medium,
                                   crowding,
                                   volume)

    particle_model.add_reaction(Reaction('r1f', {A: -1, B: -1, C: 1}, k1f))
    particle_model.add_reaction(Reaction('r1b', {A: 1, B: 1, C: -1}, k1b))

    # Define initial conditions
    particle_model.initial_conditions['A'] = parameters['A_0']
    particle_model.initial_conditions['B'] = parameters['B_0']
    particle_model.initial_conditions['C'] = parameters['C_0']

    result = particle_model.simulate(dt=parameters['dt'],
                                     max_time=parameters['t_max'],
                                     log_step=10,
                                     n_sample=0,
                                     random_seed=seed,
                                     is_hardsphere=True,
                                     is_constant_state=False,
                                     t_equlibriate=0.0)


    # Write in a data frame
    data = np.array([result.time, result.species['A'], result.species['B'], result.species['C'] ])
    df = DataFrame(data=data.T, columns = ['time', 'A', 'B', 'C'])

    return df


def ecell4_gfrd_simulation(parameters,phi= 0.0, seed=1):
    from ecell4 import Species, NetworkModel, \
        create_binding_reaction_rule, create_unbinding_reaction_rule,\
        Real3,FixedIntervalNumberObserver
    from ecell4.egfrd import EGFRDWorld, EGFRDSimulator
    from ecell4.core import GSLRandomNumberGenerator

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
        R_crw = mass2rad(parameters['mu_mass'])*1e-9
        volume_crw = 4.0 / 3.0 * np.pi * R_crw ** 3
        crw = Species('crw',R_crw*s, rad2diff(R_crw)*s2 )

    gamma = 4.0*np.pi*(parameters['r_A']+parameters['r_B'])*(parameters['D_A']+parameters['D_B'])*s3
    rescaled_keff = parameters['k_fwd']/AVOGADRO_NUMBER/1000*s3

    effective_k_binding = gamma*rescaled_keff/(gamma-rescaled_keff)

    print(rescaled_keff)
    print(gamma)
    print(effective_k_binding)
    m.add_reaction_rule(create_binding_reaction_rule(A, B, C, rescaled_keff))

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
        N_crw = round(parameters['volume']/1000*phi/volume_crw)
        w.add_molecules(crw, N_crw)

    obs = FixedIntervalNumberObserver(parameters['t_max']/1000.0, ['A', 'B', 'C'])
    sim = EGFRDSimulator(w)
    sim.run(parameters['t_max'], obs)

    # Write in a data frame
    data = np.array(obs.data())
    df = DataFrame(data=data, columns = ['time', 'A', 'B', 'C'])

    return df

def ecell4_brd_simulation(parameters,phi=0.0,seed=1):
    from ecell4 import Species, NetworkModel, \
        create_binding_reaction_rule, create_unbinding_reaction_rule,\
        Real3,FixedIntervalNumberObserver
    from ecell4.bd import BDWorld, BDSimulator
    from ecell4.core import GSLRandomNumberGenerator

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
        R_crw = mass2rad(parameters['mu_mass'])*1e-9
        volume_crw = 4.0 / 3.0 * np.pi * R_crw ** 3
        crw = Species('crw',R_crw*s, rad2diff(R_crw)*s2 )

    gamma = 4.0*np.pi*(parameters['r_A']+parameters['r_B'])*(parameters['D_A']+parameters['D_B'])*s3
    rescaled_keff = parameters['k_fwd']/AVOGADRO_NUMBER/1000*s3

    effective_k_binding = gamma*rescaled_keff/(gamma-rescaled_keff)

    print(rescaled_keff)
    print(gamma)
    print(effective_k_binding)
    m.add_reaction_rule(create_binding_reaction_rule(A, B, C, rescaled_keff))

    m.add_reaction_rule(create_unbinding_reaction_rule(C, A, B, parameters['k_fwd']*parameters['K_eq']))

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
        N_crw = round(parameters['volume']/1000*phi/volume_crw)
        print(N_crw)
        w.add_molecules(crw, N_crw)

    obs = FixedIntervalNumberObserver(parameters['t_max']/1000.0, ['A', 'B', 'C'])
    sim = BDSimulator(w)
    sim.run(parameters['t_max'], obs)

    # Write in a data frame
    data = np.array(obs.data())
    df = DataFrame(data=data, columns = ['time', 'A', 'B', 'C'])

    return df


def ecell4_ode_simulation(parameters,phi=0.0,seed=1):
    from ecell4 import Species, NetworkModel, \
        create_binding_reaction_rule, create_unbinding_reaction_rule,\
        Real3,FixedIntervalNumberObserver

    from ecell4.ode import ODEWorld, ODESimulator
    from ecell4.core import GSLRandomNumberGenerator

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

    gamma = 4.0*np.pi*(parameters['r_A']+parameters['r_B'])*(parameters['D_A']+parameters['D_B'])*s3
    rescaled_keff = parameters['k_fwd']/AVOGADRO_NUMBER/1000*s3

    effective_k_binding = gamma*rescaled_keff/(gamma-rescaled_keff)

    print(rescaled_keff)
    print(gamma)
    print(effective_k_binding)
    m.add_reaction_rule(create_binding_reaction_rule(A, B, C, rescaled_keff))

    m.add_reaction_rule(create_unbinding_reaction_rule(C, A, B, parameters['k_fwd']*parameters['K_eq']))


    a = (parameters['volume'] / 1000) ** (1 / 3) * s  # in m
    w = ODEWorld(Real3(a, a, a))
    w.bind_to(m)

    # Add the species in the concentrations
    N_A = parameters['A_0']*parameters['volume']*AVOGADRO_NUMBER
    N_B = parameters['B_0']*parameters['volume']*AVOGADRO_NUMBER
    N_C = parameters['C_0']*parameters['volume']*AVOGADRO_NUMBER


    w.add_molecules(A, N_A)
    w.add_molecules(B, N_B)
    w.add_molecules(C, N_C)


    obs = FixedIntervalNumberObserver(parameters['t_max']/1000.0, ['A', 'B', 'C'])
    sim = ODESimulator(w)
    sim.run(parameters['t_max'], obs)

    # Write in a data frame
    data = np.array(obs.data())
    df = DataFrame(data=data, columns = ['time', 'A', 'B', 'C'])

    return df



"""
Run A simulation
"""
import sys

#if __name__ is "__main__":

param_type = sys.argv[1]
sim_type   = sys.argv[2]
phi        = float(sys.argv[3])
seed       = int(sys.argv[4])
output     = sys.argv[5]


if param_type == 'diff':
    parameters = parameters_diff_lim
elif param_type == 'react':
    parameters = parameters_reaction_lim
else:
    raise ValueError('"{}" is not a valid input argument'.format(param_type))

if sim_type == 'ode':
    data = ecell4_ode_simulation(parameters,phi=phi,seed=seed)
elif sim_type == 'brd':
    data = ecell4_brd_simulation(parameters,phi=phi,seed=seed)
elif sim_type == 'gfrd':
    data = ecell4_gfrd_simulation(parameters,phi=phi,seed=seed)
elif sim_type == 'openbread':
    data = openbread_simulation(parameters,phi=phi,seed=seed)
elif sim_type == 'geek':
    data = geek_simulations(parameters,phi=phi,seed=seed)
else:
    raise ValueError('"{}" is not a valid input argument'.format(sim_type))


filename = '{}/{}_{}_{}_{}.csv'.format(output,param_type,sim_type,phi,seed)
print("Write output file {}".format(filename))

data.to_csv(filename)






