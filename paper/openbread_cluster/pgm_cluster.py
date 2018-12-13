#!$PYTHON_ROOT/bin/python3

"""
General single enzyme model to run at the cluster

"""

from skimpy.core import *
from skimpy.mechanisms import *

import sys

args = sys.argv

realization = int(args[1])
volume_fraction = float(args[2])
mu_mass = float(args[3])
sigma_mass = float(args[4])
S_concentration = float(args[5])
P_concentration = float(args[6])
E_concentration = float(args[7])
ES_concentration = float(args[8])
output = args[9]
index  = int(args[10])
job_id = int(args[11])


"""
OPENBREAD model to extract the GEEK parameters ready for use on the cluster
"""


from openbread.core import *


# Construct species in the model

g3p = Species(  name = 'g3p',
                diffusion_constant = 940e-12,
                collision_radius   = 1.11e-9,
                mass = 0. )
g2p = Species(  name = 'g2p',
                diffusion_constant = 940e-12,
                collision_radius   = 1.11e-9,
                mass = 0.186 )
pgm = Species(  name = 'pgm',
                diffusion_constant = 84.8e-12 ,
                collision_radius   = 3.87e-9,
                mass = 61 )
EC_pgm = Species(  name = 'EC_pgm',
                   diffusion_constant = 84.8e-12,
                   collision_radius   = 3.87e-9,
                   mass = 61.186 )

species = [g3p, g2p, pgm, EC_pgm]

# Define microscopic reaction rate constants:
k1f = 1e5       # 1/Ms
k1b = 20.0      # 1/s
k2f = 10.0      # 1/s
k2b = 1e5       # 1/Ms


# Setup particle simulation environemnt
volume = 10e-18 # (0.1 mum)^3 in L

medium = ParticleModel.Medium(  viscosity=0.7e-3, # Pa s
                                temperatur=310.15)

crowding = ParticleModel.Crowding( volume_fraction = volume_fraction,
                                   mu = np.log(mu_mass),
                                   sigma = sigma_mass,
                                   max_size = 10e-3)

particle_model = ParticleModel(medium,
                               crowding,
                               volume)

particle_model.add_reaction(Reaction('r1f', {g2p:-1,pgm:-1,EC_pgm:1},  k1f ))
particle_model.add_reaction(Reaction('r1b', {g2p: 1,pgm: 1,EC_pgm:-1}, k1b ))
particle_model.add_reaction(Reaction('r2f', {g3p: 1,pgm: 1,EC_pgm:-1}, k2f ))
particle_model.add_reaction(Reaction('r2b', {g3p:-1,pgm:-1,EC_pgm:1},  k2b ))


# Define initial conditions
particle_model.initial_conditions['pgm'] = 50e-6
particle_model.initial_conditions['EC_pgm'] = 50e-6
particle_model.initial_conditions['g3p'] = 50e-6
particle_model.initial_conditions['g2p'] = 50e-6



result = particle_model.simulate(   dt=0.25e-9,
                                    max_time=1e-6,
                                    log_step=10,
                                    random_seed=realization,
                                    is_hardsphere=True,
                                    is_constant_state=True,
                                    t_equlibriate=0.0)

#print("--- %s seconds ---" % (time.time() - start_time))
# Wirte the result in to a data frame and to csv

from pandas import DataFrame
columns = []
df = DataFrame(columns = columns)


map_elementray_steps = {'r1f':'k1_fwd',
                        'r1b':'k1_bwd',
                        'r2f':'k2_fwd',
                        'r2b':'k2_bwd',}

this_reaction = this_multiscale_model.reactions.pgm
this_elementrary_rate_constants = this_reaction.mechanism.elementary_reactions
this_rate_constans = this_reaction.mechanism.rate_constants
this_results = {}
this_results['realization'] = realization
this_results['mu_mass'] = mu_mass
this_results['sigma_mass'] = sigma_mass
this_results['volume_fraction'] = volume_fraction
this_results['substrate_concentration'] = S_concentration
this_results['product_concentration'] = P_concentration
this_results['enzyme_concentration'] = E_concentration
this_results['enzyme_complex_concentration'] = ES_concentration

for field in this_elementrary_rate_constants._fields:
  particle_result_key = getattr(this_elementrary_rate_constants,field)
  skimkpy_results_key = map_elementray_steps[field]
  rate_constant = result.effective_rate_constants[particle_result_key.__str__()]
  rel_rate_constant = rate_constant/getattr(this_rate_constans,skimkpy_results_key)

  this_results[skimkpy_results_key+'_effective'] = rate_constant
  this_results[skimkpy_results_key+'_relative']  = rel_rate_constant

  #key_relative_rate_constants.append(skimkpy_results_key+'_relative')


df = df.append(this_results, ignore_index=True)
df.to_csv( output  \
          +'effective_rates_var_concentration_S_P_' \
          +'phi_'+str(volume_fraction)+'_'\
          +'mu_'+str(mu_mass)+'_'\
          +'sigma_'+str(sigma_mass)\
          +'_'+str(index)\
          +'_'+str(job_id)+'.csv')
