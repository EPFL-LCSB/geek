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
    param_type,sim_type, volume_fraction, seed=re.sub('\.csv$', '', file).split('_')
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
Validation hsbrd dilute
"""

volume = 1e-18
AVOGADRO_NUMBER = 6e23

s = 1./volume/AVOGADRO_NUMBER*1e6


f, ((ax1, ax2, ), (ax3, ax4,), ) = plt.subplots(2, 2, sharey='row')


"""
Time plot
"""

phi = 0.0

for i in range(1,11):
    selection = (df['sim_type'] == 'geekcf') & \
                (df['param_type'] == 'react') & \
                (df['volume_fraction'] == phi) & \
                (df['seed'] == i)

    ax1.plot(df[selection]['time']*1e6, df[selection]['A']*s , '-', color='blue', alpha=0.5)
#    plt.plot(df[selection]['time']*1e6, df[selection]['C'], '-', color='navy', alpha=0.5)

    selection = (df['sim_type'] == 'crwdfree') & \
                (df['param_type'] == 'react') & \
                (df['volume_fraction'] == phi) & \
                (df['seed'] == i)

    ax1.plot(df[selection]['time']*1e6,df[selection]['A']*s , '-', color='firebrick', alpha=0.5)
#    plt.plot(df[selection]['time']*1e6, df[selection]['C'], '-', color='dodgerblue', alpha=0.5)



for i in range(1,11):
    selection = (df['sim_type'] == 'geekcf') & \
                (df['param_type'] == 'diff') & \
                (df['volume_fraction'] == phi) & \
                (df['seed'] == i)

    ax3.plot(df[selection]['time']*1e6, df[selection]['A']*s , '-', color='blue', alpha=0.5)
    #plt.plot(df[selection]['time']*1e6, df[selection]['C'], '-', color='navy', alpha=0.5)

    selection = (df['sim_type'] == 'crwdfree') & \
                (df['param_type'] == 'diff') & \
                (df['volume_fraction'] == phi) & \
                (df['seed'] == i)

    ax3.plot(df[selection]['time']*1e6,df[selection]['A']*s , '-', color='firebrick', alpha=0.5)
    #plt.plot(df[selection]['time']*1e6, df[selection]['C'], '-', color='dodgerblue', alpha=0.5)



#for phi in [0., 0.1, 0.2 , 0.3, 0.5]:

"""
Reaction controlled
"""

selection = (df['sim_type'] == 'crwdfree') & \
            (df['param_type'] == 'react') & \
            (df['time'] > 1e-3*0.5)

df_temp = df[selection]
df_temp['A'] = df_temp['A']*s
grouped_df = df_temp.groupby('volume_fraction')
x = grouped_df.mean().index
y = grouped_df.mean()['A'].values
yl = y - grouped_df.quantile(0.25)['A'].values
yu = grouped_df.quantile(0.75)['A'].values - y

ax2.errorbar(x, y, yerr=[yl.T,yu.T], color='firebrick', fmt='o', fillstyle='none')


selection = (df['sim_type'] == 'geekcf') & \
            (df['param_type'] == 'react') & \
            (df['time'] > 1e-3*0.5)

df_temp = df[selection]
df_temp['A'] = df_temp['A']*s
grouped_df = df_temp.groupby('volume_fraction')
x = grouped_df.mean().index
y = grouped_df.mean()['A'].values
yl = y - grouped_df.quantile(0.25)['A'].values
yu = grouped_df.quantile(0.75)['A'].values - y

ax2.errorbar(x, y, yerr=[yl.T,yu.T], color='blue', fmt='s', fillstyle='none')


"""
Difusion controlled
"""


selection = (df['sim_type'] == 'crwdfree') & \
            (df['param_type'] == 'diff') & \
            (df['time'] > 1e-5*0.5)

df_temp = df[selection]
df_temp['A'] = df_temp['A']*s
grouped_df = df_temp.groupby('volume_fraction')
x = grouped_df.mean().index
y = grouped_df.mean()['A'].values
yl = y - grouped_df.quantile(0.25)['A'].values
yu = grouped_df.quantile(0.75)['A'].values - y

ax4.errorbar(x, y, yerr=[yl.T,yu.T], color='firebrick', fmt='o', fillstyle='none')


selection = (df['sim_type'] == 'geekcf') & \
            (df['param_type'] == 'diff') & \
            (df['time'] > 1e-5*0.5)

df_temp = df[selection]
df_temp['A'] = df_temp['A']*s
grouped_df = df_temp.groupby('volume_fraction')
x = grouped_df.mean().index
y = grouped_df.mean()['A'].values
yl = y - grouped_df.quantile(0.25)['A'].values
yu = grouped_df.quantile(0.75)['A'].values - y

ax4.errorbar(x, y, yerr=[yl.T,yu.T], color='blue', fmt='s', fillstyle='none')

# Labels and legend
ax3.set_xlabel('time / [$\mu s]$')
ax3.set_ylabel('[A] / [$\mu$M]')
ax3.set_ylim([0, 55])
ax3.set_xlim([0, 10])

#plt.xlabel('time in $\mu s$')

ax1.set_ylabel('[A] / [$\mu$M]')
# plt.legend(['GEEK',
#             'hsbrd', ],
#             loc='upper center',
#             bbox_to_anchor=(0.5,1.5),
#             ncol=2)

ax1.set_ylim([0, 55])
ax1.set_xlim([0, 1000])


ax4.set_xlabel('$\phi$')

plt.tight_layout()
plt.savefig('verification_openbread_concentration.png', ppi=1200)




"""
Initial rate 
"""




def get_initial_rates(df,species):

    y = []
    dyl = []
    dyu = []

    dftemp = pd.DataFrame(columns=['v_init','volume_fraction','seed'])


    for volume_fraction in df['volume_fraction'].unique():
        for seed in df[df['volume_fraction']==volume_fraction]['seed'].unique():

            selection = (df['seed'] == seed) & \
                        (df['volume_fraction'] == volume_fraction)

            t_min = df[selection]['time'].min()*1e6
            t_max = df[selection]['time'].max()*1e6

            ix_min = df[selection]['time'].argmin()
            ix_max = df[selection]['time'].argmax()

            a0 = df[selection][species].loc[ix_min]
            a1 = df[selection][species].loc[ix_max]

            v_init = (a1 - a0) / (t_max - t_min)

            data = dict(v_init= - v_init,
                        volume_fraction=volume_fraction,
                        seed=seed)

            dftemp = dftemp.append(data, ignore_index=True)


    y = dftemp.groupby('volume_fraction').mean()['v_init'].values
    x = dftemp.groupby('volume_fraction').mean().index
    dyl = y - dftemp.groupby('volume_fraction').quantile(0.25)['v_init'].values
    dyu = dftemp.groupby('volume_fraction').quantile(0.75)['v_init'].values - y

    print(y)


    return x, y, dyl, dyu



"""Plot initial rates"""

f, (ax1, ax2, ) = plt.subplots(1,2, sharex=True, figsize=(6,3) )


selection = (df['sim_type'] == 'crwdfree') & \
            (df['param_type'] == 'react') & \
            (df['time'] < 1e-3*0.1)

x,y,dyl,dyu = get_initial_rates(df[selection], 'A')
ax1.errorbar(x, y, yerr=[dyl,dyu], color='firebrick', fmt='o', fillstyle='none')


selection = (df['sim_type'] == 'geekcf') & \
            (df['param_type'] == 'react') & \
            (df['time'] < 1e-3*0.1)

x,y,dyl,dyu = get_initial_rates(df[selection], 'A')
ax1.errorbar(x, y, yerr=[dyl,dyu], color='blue', fmt='o', fillstyle='none')




selection = (df['sim_type'] == 'crwdfree') & \
            (df['param_type'] == 'diff') & \
            (df['time'] < 1e-3*0.1)

x,y,dyl,dyu = get_initial_rates(df[selection], 'A')
ax2.errorbar(x, y, yerr=[dyl,dyu], color='firebrick', fmt='o', fillstyle='none')


selection = (df['sim_type'] == 'geekcf') & \
            (df['param_type'] == 'diff') & \
            (df['time'] < 1e-3*0.1)

x,y,dyl,dyu = get_initial_rates(df[selection], 'A')
ax2.errorbar(x, y, yerr=[dyl,dyu], color='blue', fmt='o', fillstyle='none')


ax1.set_ylabel('$\Delta$[A]/$\Delta t$ / [$\mu$M $\\mathrm{s^{-1}}$ ]')

ax1.set_xlabel('$\phi$')
ax2.set_xlabel('$\phi$')

plt.tight_layout()
plt.savefig('verification_openbread_initial_rate.png', ppi=1200)