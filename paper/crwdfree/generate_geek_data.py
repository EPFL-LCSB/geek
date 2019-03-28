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
import pandas as pd

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
    'C_0': 50e-6,       # M
    'volume': 1e-18,    # L
    't_max': 1e-6,      # s
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
    'C_0': 50e-6,       # M
    'volume': 1e-18,    # L
    't_max': 1e-6,      # s
    'dt': 1e-9,         # s
    'mu_mass': 21.1,    # s
    'du': 1e9           # kbT
}




def run_simulation(parameters,phi,seed):
    from paper.crwdfree.crowder_free_simulation import crowder_free_simulation_method
    from paper.crwdfree.crowder_free_simulation import calc_effective_volume, mass2rad, AVOGADRO_NUMBER

    dt_log = parameters['dt']*10.0

    result = crowder_free_simulation_method(parameters, phi, seed, is_geek=True, dt_log=dt_log)

    s  = float(1.0e6)
    s2 = float(1.0e12)
    s3 = float(1.0e18)

    # Reaction volume
    volume_AB = calc_effective_volume((parameters['D_A'] + parameters['D_B'])*s2,
                                      (parameters['r_A'] + parameters['r_B'])*s,
                                      dt_log)
    #
    gamma = 4.0 * np.pi * (parameters['r_A'] + parameters['r_B']) * (parameters['D_A'] + parameters['D_B']) * s3
    rescaled_keff = parameters['k_fwd'] / AVOGADRO_NUMBER / 1000.0 * s3
    effective_k_binding = gamma * rescaled_keff / (gamma - rescaled_keff)
    effective_k_unbinding = parameters['k_fwd'] * parameters['K_eq']

    n = round(0.5 * len(result.collisions))
    reaction_propablity = 1 - np.exp( -effective_k_binding*dt_log/volume_AB)

    k_1_bwd_eff_rel = np.mean(result.acceptance[n:-1])/100.0
    k_1_bwd_eff = parameters['k_fwd'] * parameters['K_eq'] * k_1_bwd_eff_rel

    k_1_fwd_eff = np.mean(result.collisions[n:-1]) / (dt_log * result.species['A'][0] * result.species['B'][0]) \
                  * reaction_propablity
    k_1_fwd_eff_rel = k_1_fwd_eff / parameters['k_fwd']

    data = np.array([[parameters['A_0'],
                     parameters['B_0'],
                     parameters['C_0'],
                     k_1_bwd_eff,
                     k_1_bwd_eff_rel,
                     k_1_fwd_eff,
                     k_1_fwd_eff_rel,
                     parameters['mu_mass'],
                     seed,
                     0.,
                     phi,
                     ],])

    df = pd.DataFrame(data.T, columns = ['A_concentration',
                                           'B_concentration',
                                           'C_concentration',
                                           'k1_bwd_effective',
                                           'k1_bwd_relative',
                                           'k1_fwd_effective',
                                           'k1_fwd_relative',
                                           'mu_mass',
                                           'realization',
                                           'sigma_mass',
                                           'volume_fraction'])
    return df



"""
Run A simulation
"""
import sys


param_type     = sys.argv[1]
phi            = float(sys.argv[2])
seed           = int(sys.argv[3])
conc_scaling_A   = float(sys.argv[4])
conc_scaling_B   = float(sys.argv[5])
conc_scaling_C   = float(sys.argv[6])
output     = sys.argv[7]


if param_type == 'diff':
    parameters = parameters_diff_lim
elif param_type == 'react':
    parameters = parameters_reaction_lim
else:
    raise ValueError('"{}" is not a valid input argument'.format(param_type))

parameters['A_0'] = parameters['A_0']*conc_scaling_A
parameters['B_0'] = parameters['B_0']*conc_scaling_B
parameters['C_0'] = parameters['C_0']*conc_scaling_C


data = run_simulation(parameters,phi=phi,seed=seed)



filename = '{}/{}_{}_{}_{}_{}_{}_{}.csv'.format(output,param_type,phi,seed,conc_scaling_A,conc_scaling_B,conc_scaling_C)
print("Write output_old file {}".format(filename))

data.to_csv(filename)


