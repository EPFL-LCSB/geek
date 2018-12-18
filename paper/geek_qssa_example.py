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
from pandas import read_csv, DataFrame
import random
import numpy as np
from geek.analysis import geek_regression


"""
Define the GEEK Michaelis-Menten ODE-System

"""

seed = 1

df = read_csv('../data/result_full_factorial_pgm.csv')

# Reference concentrations
pgm = 64e-6

total_enzyme = pgm

g3p = 49e-6
g2p = g3p

# Define microscopic reaction rate constants:
k1f = 1.52e5       # 1/Ms
k1b = 10.0      # 1/s
k2f = 22.0      # 1/s
k2b = 3.29e5       # 1/Ms



reference_concentrations = [pgm*0.5, pgm*0.5, g3p, g2p]
concentrations = ['enzyme_complex_concentration',
                  'enzyme_concentration',
                  'product_concentration',
                  'substrate_concentration']


# Filter the data frame for specific condition
this_volume_fraction = 0.1
this_mu              = 31.9
this_sigma           = 0.825

this_df = df [ (df['sigma_mass'] == this_sigma) &
               (df['mu_mass']    == this_mu) &
               (df['volume_fraction'] == this_volume_fraction)]


# Extract the GEEK parameters from Linear regression
k1_fwd_params = geek_regression(this_df,
                                  concentrations,
                                  reference_concentrations,
                                  'k1_fwd_relative',
                                  verbose=False)

k1_bwd_params = geek_regression(this_df,
                                  concentrations,
                                  reference_concentrations,
                                  'k1_bwd_relative',
                                  verbose=False)
k2_fwd_params = geek_regression(this_df,
                                  concentrations,
                                  reference_concentrations,
                                  'k1_fwd_relative',
                                  verbose=False)

k2_bwd_params = geek_regression(this_df,
                                  concentrations,
                                  reference_concentrations,
                                  'k2_bwd_relative',
                                  verbose=False)



random.seed(seed)
#Map to parameter dict
param_dict = {
    'k_1f0': k1f,
    'k_1b0': k1b,
    'beta_1f': k1_fwd_params['beta_lb'] + (k1_fwd_params['beta_ub'] - k1_fwd_params['beta_lb']) * random.random(),
    'alpha_ES_1f': k1_fwd_params['alpha_enzyme_complex_concentration_lb'] + (
            k1_fwd_params['alpha_enzyme_complex_concentration_ub'] - k1_fwd_params[
            'alpha_enzyme_complex_concentration_lb']) * random.random(),
    'alpha_E_1f': k1_fwd_params['alpha_enzyme_concentration_lb'] + (
            k1_fwd_params['alpha_enzyme_concentration_ub'] - k1_fwd_params[
        'alpha_enzyme_concentration_lb']) * random.random(),
    'alpha_P_1f': k1_fwd_params['alpha_product_concentration_lb'] + (
            k1_fwd_params['alpha_product_concentration_ub'] - k1_fwd_params[
        'alpha_product_concentration_lb']) * random.random(),
    'alpha_S_1f': k1_fwd_params['alpha_substrate_concentration_lb'] + (
            k1_fwd_params['alpha_substrate_concentration_ub'] - k1_fwd_params[
        'alpha_substrate_concentration_lb']) * random.random(),

    'beta_1b': k1_bwd_params['beta_lb'] + (k1_bwd_params['beta_ub'] - k1_bwd_params['beta_lb']) * random.random(),
    'alpha_ES_1b': k1_bwd_params['alpha_enzyme_complex_concentration_lb'] + (
            k1_bwd_params['alpha_enzyme_complex_concentration_ub'] - k1_bwd_params[
        'alpha_enzyme_complex_concentration_lb']) * random.random(),
    'alpha_E_1b': k1_bwd_params['alpha_enzyme_concentration_lb'] + (
            k1_bwd_params['alpha_enzyme_concentration_ub'] - k1_bwd_params[
        'alpha_enzyme_concentration_lb']) * random.random(),
    'alpha_P_1b': k1_bwd_params['alpha_product_concentration_lb'] + (
            k1_bwd_params['alpha_product_concentration_ub'] - k1_bwd_params[
        'alpha_product_concentration_lb']) * random.random(),
    'alpha_S_1b': k1_bwd_params['alpha_substrate_concentration_lb'] + (
            k1_bwd_params['alpha_substrate_concentration_ub'] - k1_bwd_params[
        'alpha_substrate_concentration_lb']) * random.random(),

    'k_2f0': k2f,
    'k_2b0': k2b,
    'beta_2f': k2_fwd_params['beta_lb'] + (k2_fwd_params['beta_ub'] - k2_fwd_params['beta_lb']) * random.random(),
    'alpha_ES_2f': k2_fwd_params['alpha_enzyme_complex_concentration_lb'] + (
            k2_fwd_params['alpha_enzyme_complex_concentration_ub'] - k2_fwd_params[
        'alpha_enzyme_complex_concentration_lb']) * random.random(),
    'alpha_E_2f': k2_fwd_params['alpha_enzyme_concentration_lb'] + (
            k2_fwd_params['alpha_enzyme_concentration_ub'] - k2_fwd_params[
        'alpha_enzyme_concentration_lb']) * random.random(),
    'alpha_P_2f': k2_fwd_params['alpha_product_concentration_lb'] + (
            k2_fwd_params['alpha_product_concentration_ub'] - k2_fwd_params[
        'alpha_product_concentration_lb']) * random.random(),
    'alpha_S_2f': k2_fwd_params['alpha_substrate_concentration_lb'] + (
            k2_fwd_params['alpha_substrate_concentration_ub'] - k2_fwd_params[
        'alpha_substrate_concentration_lb']) * random.random(),

    'beta_2b': k2_bwd_params['beta_lb'] + (k2_bwd_params['beta_ub'] - k2_bwd_params['beta_lb']) * random.random(),
    'alpha_ES_2b': k1_bwd_params['alpha_enzyme_complex_concentration_lb'] + (
            k2_bwd_params['alpha_enzyme_complex_concentration_ub'] - k2_bwd_params[
        'alpha_enzyme_complex_concentration_lb']) * random.random(),
    'alpha_E_2b': k2_bwd_params['alpha_enzyme_concentration_lb'] + (
            k2_bwd_params['alpha_enzyme_concentration_ub'] - k2_bwd_params[
        'alpha_enzyme_concentration_lb']) * random.random(),
    'alpha_P_2b': k2_bwd_params['alpha_product_concentration_lb'] + (
            k2_bwd_params['alpha_product_concentration_ub'] - k2_bwd_params[
        'alpha_product_concentration_lb']) * random.random(),
    'alpha_S_2b': k2_bwd_params['alpha_substrate_concentration_lb'] + (
            k2_bwd_params['alpha_substrate_concentration_ub'] - k2_bwd_params[
        'alpha_substrate_concentration_lb']) * random.random(),

    'ES0': reference_concentrations[0],
    'E0': reference_concentrations[1],
    'P0': reference_concentrations[2],
    'S0': reference_concentrations[2],
    'E_tot':total_enzyme,
}

"""
Declare ODE-Problem
"""
from sympy import symbols,Symbol, nsolve
from sympy import N as numerical
from sympy import re as real
from sympy import im as imag
from sympy import exp as sym_exp

# Variables
ES, E, P, S = symbols(['ES', 'E', 'P', 'S'])
variables = [ES, E, P, S]
# Parameters
k_1f0, k_1b0, k_2f0, k_2b0, = symbols(['k_1f0', 'k_1b0', 'k_2f0','k_2b0'] )

# Define symbols for the GEEK parameters
beta_1f, beta_1b, beta_2f, beta_2b , = symbols(['beta_1f', 'beta_1b', 'beta_2f', 'beta_2b' ] )
alpha_ES_1b,alpha_ES_1f,alpha_ES_2b,alpha_ES_2f, = symbols(['alpha_ES_1f', 'alpha_ES_1b','alpha_ES_2b','alpha_ES_2f'])
alpha_E_1b, alpha_E_1f, alpha_E_2b, alpha_E_2f, = symbols(['alpha_E_1f', 'alpha_E_1b','alpha_E_2b','alpha_E_2f'])
alpha_P_1f, alpha_P_1b, alpha_P_2f, alpha_P_2b, = symbols(['alpha_P_1f', 'alpha_P_1b','alpha_P_2f','alpha_P_2b'])
alpha_S_1f, alpha_S_1b, alpha_S_2f, alpha_S_2b, = symbols(['alpha_S_1f', 'alpha_S_1b','alpha_S_2f','alpha_S_2b'])

ES0,E0,P0, S0 = symbols(['ES0', 'E0', 'P0', 'S0'])

ode_params = [k_1f0, k_1b0, k_2f0, k_2b0,
              beta_1f, beta_1b, beta_2f, beta_2b ,
              alpha_ES_1b,alpha_ES_1f,alpha_ES_2b,alpha_ES_2f,
              alpha_E_1b, alpha_E_1f, alpha_E_2b, alpha_E_2f,
              alpha_P_1f, alpha_P_1b, alpha_P_2f, alpha_P_2b,
              alpha_S_1f, alpha_S_1b, alpha_S_2f, alpha_S_2b,
              ES0, E0, P0, S0]
# Reactions

geek_reactions = {
    'r_1f': k_1f0 * S * E * sym_exp(beta_1f) * (ES / ES0) ** alpha_ES_1f * (E / E0) ** alpha_E_1f * ( P / P0) ** alpha_P_1f * ( S / S0) ** alpha_S_1f,
    'r_1b': k_1b0 * ES * sym_exp(beta_1b) * (ES / ES0) ** alpha_ES_1b * (E / E0) ** alpha_E_1b * ( P / P0 ) ** alpha_P_1b * ( S / S0) ** alpha_S_1b,
    'r_2f': k_2f0 * ES * sym_exp(beta_2f) * (ES / ES0) ** alpha_ES_2f * (E / E0) ** alpha_E_2f * (P / P0) ** alpha_P_2f * (S / S0) ** alpha_S_2f,
    'r_2b': k_2b0 * P * E * sym_exp(beta_2b) * (ES / ES0) ** alpha_ES_2b * (E / E0) ** alpha_E_2b * (P / P0) ** alpha_P_2b * (S / S0) ** alpha_S_2b,
}

#Expressions

expressions = {
    ES: geek_reactions['r_1f']  + geek_reactions ['r_2b']- geek_reactions['r_1b'] - geek_reactions['r_2f'],
    E: -(geek_reactions['r_1f']  + geek_reactions ['r_2b']- geek_reactions['r_1b'] - geek_reactions['r_2f']),
    S: -geek_reactions['r_1f'] + geek_reactions['r_1b'],
    P: geek_reactions['r_2f'] - geek_reactions['r_2b'],
}


E_tot = symbols('E_tot')

# QSSA Expression:
ES_qss = expressions[ES].subs(E,E_tot-ES)

# Substitute values
ES_qss.subs({ Symbol(key):val for key,val in param_dict.items()})

p = np.arange(1e-2,10,0.01)*g3p
s = np.arange(1e-2,10,0.01)*g2p

data=np.array([[0., 0., 0.],])

for the_p in p:
    for the_s in s:
        ES_qss_temp = ES_qss.subs({Symbol(key): val for key, val in param_dict.items()})
        ES_qss_temp = ES_qss_temp.subs({S:the_s, P:the_p})
        ES_qss_val = nsolve(ES_qss_temp, ES, reference_concentrations[0])
        if imag(ES_qss_val) < 1e-9:
            ES_qss_val = real(ES_qss_val)
            E_qss_val = param_dict['E_tot'] - ES_qss_val
        else:
            raise ValueError('Invalid output!')


        v_net = geek_reactions['r_2f'] - geek_reactions['r_2b']

        this_params = {'ES':ES_qss_val,'E':E_qss_val,'S':the_s, 'P':the_p }
        this_params.update(param_dict)

        tamp = v_net.subs({Symbol(key): val for key, val in this_params.items()})
        v_net_val = numerical(tamp)

        data =  np.append(data,np.array([[the_s,the_p,v_net_val]]),axis=0 )

df = DataFrame(data=data, columns=['S', 'P', 'v_net'])