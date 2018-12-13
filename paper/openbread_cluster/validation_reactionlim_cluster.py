#!$PYTHON_ROOT/bin/python3

"""
General single enzyme model to run at the cluster

"""
import sys

args = sys.argv

realization = int(args[1])
volume_fraction = float(args[2])
mu_mass = float(args[3])
sigma_mass = float(args[4])
A_concentration = float(args[5])
B_concentration = float(args[6])
C_concentration = float(args[7])

output = args[8]
index  = int(args[9])
job_id = int(args[10])


"""
OPENBREAD model to extract the GEEK parameters ready for use on the cluster
"""

# parameters_reaction_lim = {
#     'K_eq':  50e-6,     # M
#     'k_fwd': 5e7,       # 1/Ms
#     'r_A': 2e-9,        # m
#     'r_B': 2e-9,        # m
#     'r_C': 3e-9,        # m
#     'D_A': 500e-12,     # m^2/s
#     'D_B': 500e-12,     # m^2/s
#     'D_C': 350e-12,     # m^2/s
#     'm_A': 10,          # kDa
#     'm_B': 10,          # kDa
#     'm_C': 20,          # kDa
#     'A_0': 50e-6,       # M
#     'B_0': 50e-6,       # M
#     'C_0': 0,           # M
#     'volume': 10e-18,   # L
#     't_max': 1e-3,      # s
#     'dt': 0.25e-9,      # s
#     'mu_mass': 21.1     # s
# }



from openbread.core import *


# Construct species in the model

A = Species(name = 'A',
            diffusion_constant = 500e-12,
            collision_radius   = 2e-9,
            mass = 10 )
B = Species(name = 'B',
            diffusion_constant = 500e-12,
            collision_radius   = 2e-9,
            mass = 10 )

C = Species(name = 'C',
            diffusion_constant = 350e-12,
            collision_radius   = 3e-9,
            mass = 20 )


species = [A, B, C, ]

# Define microscopic reaction rate constants:
k1f = 5e7         # 1/Ms
k1b = 5e7*50e-6   # 1/s


# Setup particle simulation environemnt
volume = 10e-18 # (0.1 mum)^3 in L

medium = ParticleModel.Medium( viscosity=0.7e-3, # Pa s
                               temperatur=310.15)

crowding = ParticleModel.Crowding( volume_fraction = volume_fraction,
                                   mu = np.log(mu_mass),
                                   sigma = sigma_mass,
                                   max_size = 10e-3)

particle_model = ParticleModel(medium,
                               crowding,
                               volume)

particle_model.add_reaction(Reaction('r1f', {A:-1, B: -1, C:1},  k1f))
particle_model.add_reaction(Reaction('r1b', {A: 1, B:  1, C:-1}, k1b))


# Define initial conditions
particle_model.initial_conditions['A'] = A_concentration
particle_model.initial_conditions['B'] = B_concentration
particle_model.initial_conditions['C'] = C_concentration


result = particle_model.simulate(   dt=0.25e-9,
                                    max_time=1e-6,
                                    log_step=10,
                                    n_sample=100,
                                    random_seed=realization,
                                    is_hardsphere=True,
                                    is_constant_state=True,
                                    t_equlibriate=0.0)

#print("--- %s seconds ---" % (time.time() - start_time))
# Wirte the result in to a data frame and to csv

from pandas import DataFrame
columns = []
df = DataFrame(columns=columns)


map_elementray_steps = {'r1f': 'k1_fwd',
                        'r1b': 'k1_bwd'
                        }


this_results = dict()
this_results['realization'] = realization
this_results['mu_mass'] = mu_mass
this_results['sigma_mass'] = sigma_mass
this_results['volume_fraction'] = volume_fraction
this_results['A_concentration'] = A_concentration
this_results['B_concentration'] = B_concentration
this_results['C_concentration'] = C_concentration

for this_reaction in map_elementray_steps.keys():

    elementary_results_key = map_elementray_steps[this_reaction]

    rate_constant = result.effective_rate_constants[this_reaction]

    rel_rate_constant = rate_constant/particle_model.reactions[this_reaction].rate_constant

    this_results[elementary_results_key + '_effective'] = rate_constant
    this_results[elementary_results_key + '_relative'] = rel_rate_constant


df = df.append(this_results, ignore_index=True)
df.to_csv(  output
            + '/effective_rates_var_concentration_A_B_'
            + 'phi_'+str(volume_fraction)+'_'
            + 'mu_'+str(mu_mass)+'_'
            + 'sigma_'+str(sigma_mass)
            + '_'+str(index)
            + '_'+str(job_id)+'.csv')
