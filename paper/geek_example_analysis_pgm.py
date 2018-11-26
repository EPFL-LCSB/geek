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

from pandas import read_csv
from geek.analysis import geek_regression

df = read_csv('../data/result_full_factorial_pgm.csv')

# Reference concentrations
pgm = 64e-6
g3p = 49e-6
g2p = g3p

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


geek_parameters = geek_regression(this_df,
                                  concentrations,
                                  reference_concentrations,
                                  'k1_bwd_relative',
                                  verbose=True)