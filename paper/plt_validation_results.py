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

import pandas as pd
import re

from os import listdir

from sys import argv



def find_csv_filenames( path_to_dir, suffix=".csv" ):
    filenames = listdir(path_to_dir)
    return [ filename for filename in filenames if filename.endswith( suffix ) ]

# Find all csv files in input path
files = find_csv_filenames(argv[1])

df = pd.DataFrame(columns=['sim_type', 'param_type', 'volume_fraction', 'seed', 'time', 'A', 'B', 'C'])

for file in files:
    param_type,sim_type,volume_fraction,seed =re.sub('\.csv$','',file).split('_')
    print(param_type)
    this_df =  pd.read_csv(argv[1]+'/'+file, index_col=0)
    this_df['sim_type'] = sim_type
    this_df['param_type'] = param_type
    this_df['volume_fraction'] = float(volume_fraction)
    this_df['seed'] = int(seed)
    df = df.append(this_df, ignore_index=True)



import matplotlib as mpl
mpl.use("pgf")
pgf_with_rc_fonts = {
    "font.family": "serif",
    "font.size": 10,
    "font.serif": [],                   # use latex default serif font
}
mpl.rcParams.update(pgf_with_rc_fonts)

import matplotlib.pyplot as plt


"""
Validation openbread dilute
"""



f = plt.figure()
phi = 0.0
for i in range(1,11):
    selection = (df['sim_type'] == 'openbread') & \
                (df['param_type'] == 'diff') & \
                (df['volume_fraction'] == phi) & \
                (df['seed'] == i)

    plt.plot(df[selection]['time']*1e6,df[selection]['A'] , '-', color='firebrick', alpha=0.5)
    plt.plot(df[selection]['time']*1e6, df[selection]['C'], '-', color='dodgerblue', alpha=0.5)

    selection = (df['sim_type'] == 'brd') & \
                (df['param_type'] == 'diff') & \
                (df['volume_fraction'] == phi) & \
                (df['seed'] == i)

    plt.plot(df[selection]['time']*1e6, df[selection]['A'] , '-', color='grey', alpha=0.5)
    plt.plot(df[selection]['time']*1e6, df[selection]['C'], '-', color='darkgrey', alpha=0.5)

# for i in range(1,6):
#     selection = (this_df['sim_type'] == 'gfrd') & \
#                 (this_df['param_type'] == 'diff') & \
#                 (this_df['volume_fraction'] == phi) & \
#                 (this_df['seed'] == i)
#
#     plt.plot(df[selection][time],df[selection]['A'] , '-', color='dogerblue', alpha=0.1)

selection = (df['sim_type'] == 'ode') & \
            (df['param_type'] == 'diff') & \
            (df['volume_fraction'] == phi) & \
            (df['seed'] == 1)

plt.plot(df[selection]['time']*1e6,df[selection]['A'] , '--', color='black')
plt.plot(df[selection]['time']*1e6,df[selection]['C'] , ':', color='black')

plt.xlabel('time in $\mu s$')
plt.ylabel('Number of molecules')
plt.legend(['A openbread', 'C openbread','A RBD', 'C RBD','A Massaction', 'C Massaction'  ])

plt.savefig('validation_openbread_dilute.png')

