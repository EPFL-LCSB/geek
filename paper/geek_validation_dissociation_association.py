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

from ecell4.bd import BDWorld, BDSimulator
from ecell4.egfrd import EGFRDWorld, EGFRDSimulator
from ecell4 import species_attributes, reaction_rules, get_model

"""
The following script shall provide a validation of GEEK and the brownian reaction dynamics by comparing its results with 
other approaches
for uni-molecular and for bi-molecular reaction
"""

"""
Diffusion limited conditions
"""

# gamma = 4*pi*4e-9*200e-12*1000*6e23 ~ 6e9

parameters = {
    'Keq':   0.8,
    'k_fwd': 1e10,      # 1/Ms
    'r_A': 2e-9,        # m
    'r_B': 2e-9,        # m
    'r_C': 3e-9,        # m
    'D_A': 100e-12,     # m^2/s
    'D_B': 100e-12,     # m^2/s
    'D_C':  75e-12,     # m^2/s
    'm_A': 10,          # kDa
    'm_B': 10,          # kDa
    'm_C': 20,          # kDa
    'A_0': 1e-6,        # M
    'B_0': 1e-6,        # M
    'volume': 10e-18,   # L
    't_max': 1e-6,      # s
    'dt': 0.25e-9,      # s
}




"""
Reaction limited conditions
"""

# gamma = 4*pi*4e-9*200e-12*1000*6e23 ~ 6e9

parameters = {
    'Keq':   0.8,
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
    'A_0': 1e-6,        # M
    'B_0': 1e-6,        # M
    'volume': 10e-18,   # L
    't_max': 1e-6,      # s
    'dt': 0.25e-9,      # s
}



def geek_simulations(parameters,phi= 0.0):
    pass


def openbread_simulation(parameters,phi= 0.0):
    pass


def ecell4_gfrd_simulation(parameters,phi= 0.0):

    with species_attributes():
        A | {'D': parameters['D_A'],  'radius':  parameters['r_A']}
        B | {'D': parameters['D_B'],  'radius':  parameters['r_B']}
        C | {'D': parameters['D_C'],  'radius':  parameters['r_C']}

    with reaction_rules():
        A + B == C | (parameters['k_fwd'],parameters['k_fwd']*parameters['K_eq'])

    m = get_model()

    w= EGFRDWorld()
    w.bind_to(m)

    obs = FixedIntervalNumberObserver(parameters['dt'], ['A', 'B', 'C'])
    sim = EGFRDSimulator(w)

    sim.run(parameters['t_max'], obs)
    return obs.data


def ecell4_brd_simulation(parameters,phi=0.0):

    with species_attributes():
        A | {'D': parameters['D_A'], 'radius': parameters['r_A']}
        B | {'D': parameters['D_B'], 'radius': parameters['r_B']}
        C | {'D': parameters['D_C'], 'radius': parameters['r_C']}

    with reaction_rules():
        A + B == C | (parameters['k_fwd'], parameters['k_fwd'] * parameters['K_eq'])

    m = get_model()
    w = BDWorld()
    w.bind_to(m)

    obs = FixedIntervalNumberObserver(parameters['dt'], ['A', 'B', 'C'])
    sim = BDSimulator(w)

    set_dt(parameters['dt'])

    sim.run(parameters['t_max'], obs)
    return obs.data





