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

from sympy import symbols
from scipy.integrate import odes

from geek.analysis import OdeFun





"""
Define the ideal Michaelis-Menten ODE-System
"""

# Variables
E, ES, S, P = symbols('E', 'ES', 'S', 'P')
variables = [E, ES, S, P]

# Parameters
k_1f0, k_1b0, k_2f0, k_2b0 = symbols('k_1f0', 'k_1b0', 'k_2f0', 'k_2b0')

parameters = [k_1f0, k_1b0, k_2f0, k_2b0]

# Reactions
reactions = {
    'r_1f': k_1f0*E*S,
    'r_1b': k_1b0*ES,
    'r_2f': k_2f0*ES,
    'r_2b': k_2b0*E*P,
}

# Expressions
expressions = {
    S:  -reactions['r_1f'] + reactions['r_1b'],
    E:  -reactions['r_1f'] + reactions['r_1b'] + reactions['r_2f'] - reactions['r_2b'],
    ES:  reactions['r_1f'] - reactions['r_1b'] - reactions['r_2f'] + reactions['r_2b'],
    P:   reactions['r_2f'] - reactions['r_2b']
}


#Add mass conservation of the enzyme

Etot = symbols('Etot')

dESdt = expressions[ES].subs(E, Etot-ES)



"""
Define the GEEK Michaelis-Menten ODE-System

"""

# Define symbols for the GEEK parameters
beta_1f, beta_1b, beta_2f, beta_2b = symbols('beta_1f', 'beta_1b', 'beta_2f', 'beta_2b')
alpha_E_1f, alpha_E_1b, alpha_E_2f, alpha_E_2b = symbols('alpha_E_1f', 'alpha_E_1b', 'alpha_E_2f', 'alpha_E_2b')
alpha_ES_1f, alpha_ES_1b, alpha_ES_2f, alpha_ES_2b = symbols('alpha_ES_1f', 'alpha_ES_1b', 'alpha_ES_2f', 'alpha_ES_2b')
alpha_S_1f, alpha_S_1b, alpha_S_2f, alpha_S_2b = symbols('alpha_S_1f', 'alpha_S_1b', 'alpha_S_2f', 'alpha_S_2b')
alpha_P_1f, alpha_P_1b, alpha_P_2f, alpha_P_2b = symbols('alpha_P_1f', 'alpha_P_1b', 'alpha_P_2f', 'alpha_P_2b')


geek_reactions = {
    'r_1f': k_1f0*E*S *exp(beta_1f)*(E/E0)**alpha_E_1f*(ES/ES0)**alpha_ES_1f*(S/S0)**alpha_P_1f*(P/P0)**alpha_P_1f,
    'r_1b': k_1b0*ES  *exp(beta_1b)*(E/E0)**alpha_E_1b*(ES/ES0)**alpha_ES_1b*(S/S0)**alpha_P_1b*(P/P0)**alpha_P_1b,
    'r_2f': k_2f0*ES  *exp(beta_2f)*(E/E0)**alpha_E_2f*(ES/ES0)**alpha_ES_2f*(S/S0)**alpha_P_2f*(P/P0)**alpha_P_2f,
    'r_2b': k_2b0*E*P *exp(beta_2b)*(E/E0)**alpha_E_2b*(ES/ES0)**alpha_ES_2b*(S/S0)**alpha_P_2b*(P/P0)**alpha_P_2b,
}

# Expressions
expressions = {
    S:  -geek_reactions['r_1f'] + geek_reactions['r_1b'],
    E:  -geek_reactions['r_1f'] + geek_reactions['r_1b'] + geek_reactions['r_2f'] - geek_reactions['r_2b'],
    ES:  geek_reactions['r_1f'] - geek_reactions['r_1b'] - geek_reactions['r_2f'] + geek_reactions['r_2b'],
    P:   geek_reactions['r_2f'] - geek_reactions['r_2b']
}

Etot = symbols('Etot')

dESdt_geek = expressions[ES].subs(E, Etot-ES)


# Get geek parameters
subs_dict = {

}