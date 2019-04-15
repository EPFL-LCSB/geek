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

import time as tim
import numpy as np
import random

from pandas import DataFrame, read_csv



AVOGADRO_NUMBER = 6e23

"""
The following script shall provide a validation of GEEK and the brownian reaction dynamics by comparing its results with 
other approaches for uni-molecular and for bi-molecular reaction

We show that the geek frame work is able to capture the behaviour by different simulation techniques



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
    'volume': 1e-18,    # L
    't_max': 1e-5,      # s
    'dt': 1.e-9,        # s
    'mu_mass': 21.1,    # s
    'du': 1e9           # kbT
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
    'volume': 1e-18,    # L
    't_max': 1e-3,      # s
    'dt': 1e-9,         # s
    'mu_mass': 21.1,    # s
    'du': 1e9           # kbT
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

from numpy import pi,sqrt,exp
from scipy.special import erfc

def calc_effective_volume(diffusion, dist, delta_t):
    """ Normalization factor for Brownian dyanamics simulation """

    # Bi mol rxn scaling
    sig_2 = 4.0 * delta_t * diffusion
    sig = sqrt(sig_2)

    exp_4_r_sig = exp(-4.0 * dist ** 2 / sig_2)

    # Expresion
    A = (sig ** 3 - 2.0 * dist ** 2 * sig) * exp_4_r_sig
    B = 6.0 * dist ** 2 * sig - sig ** 3 + 4.0 * sqrt(pi) * dist ** 3 * erfc(2.0 * dist / sig)

    effective_volume = 4.0 * pi * (A + B) / 12.0 / sqrt(pi)

    return effective_volume


""" 
Simulation functions
"""

def geek_simulations_hardsphere(parameters, sim_type, phi= 0.0, seed=1):
    from geek.analysis import geek_regression

    if sim_type == 'diff':
        df = read_csv('../data/validation_diffusion_lim_hardsphere.csv')
    elif sim_type == 'react':
        df = read_csv('../data/validation_reaction_lim_hardsphere.csv')
    else:
        raise ValueError('{} is not a valid input'.format(sim_type))
    # Reference concentrations
    reference_concentrations = [50e-6,]*3
    concentrations = ['A_concentration',
                      'B_concentration',
                      'C_concentration',]


    this_df = df[(df['volume_fraction'] == phi)]

    # Extract the GEEK parameters from Linear regression
    k1_fwd_params = geek_regression(this_df,
                                      concentrations,
                                      reference_concentrations,
                                      'k1_fwd_relative',
                                      verbose=True)

    k1_bwd_params = geek_regression(this_df,
                                      concentrations,
                                      reference_concentrations,
                                      'k1_bwd_relative',
                                      verbose=True)

    random.seed(seed)
    #Map to parameter dict
    param_dict = {
        'k_1f0': parameters['k_fwd'],
        'k_1b0': parameters['k_fwd']*parameters['K_eq'],
        'beta_1f': k1_fwd_params['beta_lb'] + (k1_fwd_params['beta_ub'] - k1_fwd_params['beta_lb']) * random.random(),
        'alpha_A_1f': k1_fwd_params['alpha_A_concentration_lb'] + (
                k1_fwd_params['alpha_A_concentration_ub'] - k1_fwd_params[
                'alpha_A_concentration_lb']) * random.random(),
        'alpha_B_1f': k1_fwd_params['alpha_B_concentration_lb'] + (
                k1_fwd_params['alpha_B_concentration_ub'] - k1_fwd_params[
            'alpha_B_concentration_lb']) * random.random(),
        'alpha_C_1f': k1_fwd_params['alpha_C_concentration_lb'] + (
                k1_fwd_params['alpha_C_concentration_ub'] - k1_fwd_params[
            'alpha_C_concentration_lb']) * random.random(),
        'beta_1b': k1_bwd_params['beta_lb'] + (k1_bwd_params['beta_ub'] - k1_bwd_params['beta_lb']) * random.random(),
        'alpha_A_1b': k1_bwd_params['alpha_A_concentration_lb'] + (
                k1_bwd_params['alpha_A_concentration_ub'] - k1_bwd_params[
            'alpha_A_concentration_lb']) * random.random(),
        'alpha_B_1b': k1_bwd_params['alpha_B_concentration_lb'] + (
                k1_bwd_params['alpha_B_concentration_ub'] - k1_bwd_params[
            'alpha_B_concentration_lb']) * random.random(),
        'alpha_C_1b': k1_bwd_params['alpha_C_concentration_lb'] + (
                k1_bwd_params['alpha_C_concentration_ub'] - k1_bwd_params[
            'alpha_C_concentration_lb']) * random.random(),
        'A0': reference_concentrations[0],
        'B0': reference_concentrations[1],
        'C0': reference_concentrations[2],
    }

    """
    Declare ODE-Problem
    """
    from sympy import symbols
    from sympy import exp as sym_exp

    # Variables
    A, B, C = symbols(['A', 'B', 'C'])
    variables = [A, B, C,]
    # Parameters
    k_1f0, k_1b0, = symbols(['k_1f0', 'k_1b0',] )
    # Define symbols for the GEEK parameters
    beta_1f, beta_1b,  = symbols(['beta_1f', 'beta_1b',] )
    alpha_A_1f, alpha_A_1b, = symbols(['alpha_A_1f', 'alpha_A_1b',])
    alpha_B_1f, alpha_B_1b, = symbols(['alpha_B_1f', 'alpha_B_1b',])
    alpha_C_1f, alpha_C_1b, = symbols(['alpha_C_1f', 'alpha_C_1b',])
    A0,B0,C0 = symbols(['A0', 'B0', 'C0'])

    ode_params = [k_1f0, k_1b0,
                  beta_1f, beta_1b ,
                  alpha_A_1b, alpha_A_1f ,
                  alpha_B_1b, alpha_B_1f,
                  alpha_C_1f, alpha_C_1b,
                  A0, B0, C0]
    # Reactions

    geek_reactions = {
        'r_1f': k_1f0 * A * B * sym_exp(beta_1f) * (A / A0) ** alpha_A_1f * (B / B0) ** alpha_B_1f * (
                    C / C0) ** alpha_C_1f,
        'r_1b': k_1b0 * C * sym_exp(beta_1b) * (A / A0) ** alpha_A_1b * (B / B0) ** alpha_B_1b * (
                    C / C0) ** alpha_C_1b
    }

    #Expressions

    expressions = {
        A: geek_reactions['r_1b'] - geek_reactions['r_1f'],
        B: geek_reactions['r_1b'] - geek_reactions['r_1f'],
        C: geek_reactions['r_1f'] - geek_reactions['r_1b'],
    }

    from geek.analysis.ode_function import OdeFun
    fun = OdeFun(variables,ode_params,expressions)

    from scipy.integrate import ode
    r = ode(fun).set_integrator('vode', method='bdf')

    eps = 1e-3

    A0 = round(parameters['A_0']*AVOGADRO_NUMBER*parameters['volume'])/AVOGADRO_NUMBER/parameters['volume']
    B0 = round(parameters['B_0']*AVOGADRO_NUMBER*parameters['volume'])/AVOGADRO_NUMBER/parameters['volume']
    C0 = round(parameters['C_0']*AVOGADRO_NUMBER*parameters['volume'])/AVOGADRO_NUMBER/parameters['volume']

    print(A0,B0,C0)

    y0 = [A0 * (1. - eps),
          B0 * (1. - eps),
          A0 * eps]
    t0 = 0.0

    r.set_initial_value(y0, t0).set_f_params(param_dict)
    data = []

    scale = parameters['volume']*AVOGADRO_NUMBER
    while r.successful() and r.t < parameters['t_max']:
        data.append( np.append(r.t + parameters['t_max']/1000.0,
                     r.integrate(r.t + parameters['t_max']/1000.0) * scale))
    data = np.array(data)

    df = DataFrame(data=data, columns = ['time', 'A', 'B', 'C'])
    return df

 
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
                                      max_size=3e-3) # For this model the max size is 3nm
    
    
    particle_model = ParticleModel(medium,
                                   crowding,
                                   volume)

    particle_model.add_reaction(Reaction('r1f', {A: -1, B: -1, C: 1}, k1f))
    particle_model.add_reaction(Reaction('r1b', {A: 1, B: 1, C: -1}, k1b))

    # Define initial conditions
    particle_model.initial_conditions['A'] = parameters['A_0']
    particle_model.initial_conditions['B'] = parameters['B_0']
    particle_model.initial_conditions['C'] = parameters['C_0']

    Nt = parameters['t_max']/parameters['dt']


    result = particle_model.simulate(dt=parameters['dt'],
                                     max_time=parameters['t_max'],
                                     log_step=max(int(Nt/1000),1),
                                     n_sample=0,
                                     random_seed=seed,
                                     is_hardsphere=True,
                                     is_constant_state=False,
                                     t_equlibriate=0.0)


    # Write in a data frame
    data = np.array([result.time, result.species['A'], result.species['B'], result.species['C'] ])
    df = DataFrame(data=data.T, columns = ['time', 'A', 'B', 'C'])

    return df


def crowder_free_simulation(parameters, phi=0.0, seed=1):
    from paper.crwdfree.crowder_free_simulation import crowder_free_simulation_method, particle, check_collision

    result = crowder_free_simulation_method(parameters, phi, seed)

    # Write in a data frame
    data = np.array([result.time, result.species['A'], result.species['B'], result.species['C'] ])
    df = DataFrame(data=data.T, columns = ['time', 'A', 'B', 'C'])

    return df


def geek_simulations_crwderfree(parameters, sim_type, phi= 0.0, seed=1):
    from geek.analysis import geek_regression

    if sim_type == 'diff':
        df = read_csv('../data/validation_diffusion_lim_crowderfree.csv')
    elif sim_type == 'react':
        df = read_csv('../data/validation_reaction_lim_crowderfree.csv')
    else:
        raise ValueError('{} is not a valid input'.format(sim_type))
    # Reference concentrations
    reference_concentrations = [50e-6,]*3
    concentrations = ['A_concentration',
                      'B_concentration',
                      'C_concentration',]


    this_df = df[(df['volume_fraction'] == phi)]

    # Extract the GEEK parameters from Linear regression
    k1_fwd_params = geek_regression(this_df,
                                      concentrations,
                                      reference_concentrations,
                                      'k1_fwd_relative',
                                      verbose=True)

    k1_bwd_params = geek_regression(this_df,
                                      concentrations,
                                      reference_concentrations,
                                      'k1_bwd_relative',
                                      verbose=True)

    random.seed(seed)
    #Map to parameter dict
    param_dict = {
        'k_1f0': parameters['k_fwd'],
        'k_1b0': parameters['k_fwd']*parameters['K_eq'],
        'beta_1f': k1_fwd_params['beta_lb'] + (k1_fwd_params['beta_ub'] - k1_fwd_params['beta_lb']) * random.random(),
        'alpha_A_1f': k1_fwd_params['alpha_A_concentration_lb'] + (
                k1_fwd_params['alpha_A_concentration_ub'] - k1_fwd_params[
                'alpha_A_concentration_lb']) * random.random(),
        'alpha_B_1f': k1_fwd_params['alpha_B_concentration_lb'] + (
                k1_fwd_params['alpha_B_concentration_ub'] - k1_fwd_params[
            'alpha_B_concentration_lb']) * random.random(),
        'alpha_C_1f': k1_fwd_params['alpha_C_concentration_lb'] + (
                k1_fwd_params['alpha_C_concentration_ub'] - k1_fwd_params[
            'alpha_C_concentration_lb']) * random.random(),
        'beta_1b': k1_bwd_params['beta_lb'] + (k1_bwd_params['beta_ub'] - k1_bwd_params['beta_lb']) * random.random(),
        'alpha_A_1b': k1_bwd_params['alpha_A_concentration_lb'] + (
                k1_bwd_params['alpha_A_concentration_ub'] - k1_bwd_params[
            'alpha_A_concentration_lb']) * random.random(),
        'alpha_B_1b': k1_bwd_params['alpha_B_concentration_lb'] + (
                k1_bwd_params['alpha_B_concentration_ub'] - k1_bwd_params[
            'alpha_B_concentration_lb']) * random.random(),
        'alpha_C_1b': k1_bwd_params['alpha_C_concentration_lb'] + (
                k1_bwd_params['alpha_C_concentration_ub'] - k1_bwd_params[
            'alpha_C_concentration_lb']) * random.random(),
        'A0': reference_concentrations[0],
        'B0': reference_concentrations[1],
        'C0': reference_concentrations[2],
    }

    """
    Declare ODE-Problem
    """
    from sympy import symbols
    from sympy import exp as sym_exp

    # Variables
    A, B, C = symbols(['A', 'B', 'C'])
    variables = [A, B, C,]
    # Parameters
    k_1f0, k_1b0, = symbols(['k_1f0', 'k_1b0',] )
    # Define symbols for the GEEK parameters
    beta_1f, beta_1b,  = symbols(['beta_1f', 'beta_1b',] )
    alpha_A_1f, alpha_A_1b, = symbols(['alpha_A_1f', 'alpha_A_1b',])
    alpha_B_1f, alpha_B_1b, = symbols(['alpha_B_1f', 'alpha_B_1b',])
    alpha_C_1f, alpha_C_1b, = symbols(['alpha_C_1f', 'alpha_C_1b',])
    A0,B0,C0 = symbols(['A0', 'B0', 'C0'])

    ode_params = [k_1f0, k_1b0,
                  beta_1f, beta_1b ,
                  alpha_A_1b, alpha_A_1f ,
                  alpha_B_1b, alpha_B_1f,
                  alpha_C_1f, alpha_C_1b,
                  A0, B0, C0]
    # Reactions

    geek_reactions = {
        'r_1f': k_1f0 * A * B * sym_exp(beta_1f) * (A / A0) ** alpha_A_1f * (B / B0) ** alpha_B_1f * (
                    C / C0) ** alpha_C_1f,
        'r_1b': k_1b0 * C * sym_exp(beta_1b) * (A / A0) ** alpha_A_1b * (B / B0) ** alpha_B_1b * (
                    C / C0) ** alpha_C_1b
    }

    #Expressions

    expressions = {
        A: geek_reactions['r_1b'] - geek_reactions['r_1f'],
        B: geek_reactions['r_1b'] - geek_reactions['r_1f'],
        C: geek_reactions['r_1f'] - geek_reactions['r_1b'],
    }

    from geek.analysis.ode_function import OdeFun
    fun = OdeFun(variables,ode_params,expressions)

    from scipy.integrate import ode
    r = ode(fun).set_integrator('vode', method='bdf')

    eps = 1e-3

    A0 = round(parameters['A_0']*AVOGADRO_NUMBER*parameters['volume'])/AVOGADRO_NUMBER/parameters['volume']
    B0 = round(parameters['B_0']*AVOGADRO_NUMBER*parameters['volume'])/AVOGADRO_NUMBER/parameters['volume']
    C0 = round(parameters['C_0']*AVOGADRO_NUMBER*parameters['volume'])/AVOGADRO_NUMBER/parameters['volume']

    print(A0,B0,C0)

    y0 = [A0 * (1. - eps),
          B0 * (1. - eps),
          A0 * eps]
    t0 = 0.0

    r.set_initial_value(y0, t0).set_f_params(param_dict)
    data = []

    scale = parameters['volume']*AVOGADRO_NUMBER
    while r.successful() and r.t < parameters['t_max']:
        data.append( np.append(r.t + parameters['t_max']/1000.0,
                     r.integrate(r.t + parameters['t_max']/1000.0) * scale))
    data = np.array(data)

    df = DataFrame(data=data, columns = ['time', 'A', 'B', 'C'])
    return df


"""
Run A simulation
"""
import sys


#if __name__ is "__main__":

param_type     = sys.argv[1]
sim_type       = sys.argv[2]
phi            = float(sys.argv[3])
seed           = int(sys.argv[4])
time_scaling   = float(sys.argv[5])
conc_scaling   = float(sys.argv[6])
dt_scaling   = float(sys.argv[7])
output     = sys.argv[8]


if param_type == 'diff':
    parameters = parameters_diff_lim
elif param_type == 'react':
    parameters = parameters_reaction_lim
else:
    raise ValueError('"{}" is not a valid input argument'.format(param_type))

parameters['t_max'] = parameters['t_max']*time_scaling
parameters['A_0'] = parameters['A_0']*conc_scaling
parameters['B_0'] = parameters['B_0']*conc_scaling
parameters['dt'] = parameters['dt']*dt_scaling

if sim_type == 'hsbrd':
    data = openbread_simulation(parameters,phi=phi,seed=seed)
elif sim_type == 'geekhs':
    data = geek_simulations_hardsphere(parameters, param_type, phi=phi, seed=seed)
elif sim_type == 'crwdfree':
    data = crowder_free_simulation(parameters, phi=phi, seed=seed)
elif sim_type == 'geekcf':
    data = geek_simulations_crwderfree(parameters, param_type, phi=phi, seed=seed)

else:
    raise ValueError('"{}" is not a valid input argument'.format(sim_type))

if time_scaling == 1 and conc_scaling == 1 and dt_scaling == 1:
    filename = '{}/{}_{}_{}_{}.csv'.format(output, param_type, sim_type, phi, seed)
else:
    filename = '{}/{}_{}_{}_{}_{}_{}_{}.csv'.format(output,param_type,sim_type,phi,seed,time_scaling,conc_scaling,dt_scaling)
print("Write output_old file {}".format(filename))

data.to_csv(filename)







