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
import statsmodels.api as statsmodels
from pandas import  DataFrame

def geek_regression(df,
                    concentrations,
                    reference_concentrations,
                    rate_constant,
                    verbose=True,
                    confidence_level=0.05):
    """

    :param df: dataframe with measured rate constant data
    :param concentrations: list of column names of the concentrations
    :param reference_concentrations: numpy array of reference concentrations
    :param rate_constant: column name of the rate constant data
    :param verbose: True or False print regression output_old
    :return: A dict with the geek parameters
    """

    if len(concentrations) != len(reference_concentrations):
        raise ValueError("You must provide a reference concentration for each concentration column")

    # Prepare multivariate linear regression
    X = np.log(df[concentrations]/reference_concentrations)
    Y = np.log(df[rate_constant])
    X = statsmodels.add_constant(X)

    # Get the odinary lineast squares residuals
    ols_result = statsmodels.OLS(Y, X).fit()
    ols_resid = DataFrame(dict(fittedvalues=ols_result.fittedvalues,
                               residuals=ols_result.resid))

    # Group the values by their fitted value -> split in 10 groups
    ols_resid['value_groups'] = (ols_resid['fittedvalues'] / max(abs(ols_resid['fittedvalues']))).round(1)

    num = 9.0
    #Count groups and check weather minimum 3 entries are in one group
    while min(ols_resid.groupby('value_groups')['residuals'].size()) < 3:
        ols_resid['value_groups'] = (ols_resid['fittedvalues'] / max(abs(ols_resid['fittedvalues'])) * num ).round(0) / num
        num -= 1.0
        if num < 1:
            raise ValueError("To little data")

    # Calculate the conditional standard deviation of the fitted values
    conditional_std = ols_resid.groupby('value_groups')['residuals'].std()
    ols_resid['std'] = np.repeat(np.NaN, len(ols_resid))
    for this_group in conditional_std.index.values:
        ols_resid['std'][ols_resid['value_groups'] == this_group] = conditional_std.ix[this_group]

    # Fit with weights antiproportinal to the standart deviation
    weights = 1.0/ols_resid['std']

    # Fit the model
    this_lin_model = statsmodels.WLS(Y, X, weights=weights)
    this_result = this_lin_model.fit()

    if verbose:
        print(this_result.summary())

    # Data for the
    if this_result.pvalues['const'] < confidence_level:
        data = {'rate_constant': rate_constant,
                'beta': this_result.params['const'],
                'beta_lb': this_result.conf_int(0.05)[0]['const'],
                'beta_ub': this_result.conf_int(0.05)[1]['const'],
                'beta_p': this_result.pvalues['const'],
                }
    else:
        data = {'rate_constant': rate_constant,
                'beta': .0,
                'beta_lb': .0,
                'beta_ub': .0,
                'beta_p': this_result.pvalues['const'],
                }

    for this_conc in concentrations:
        if this_result.pvalues[this_conc] < confidence_level:
            data['alpha_{}'.format(this_conc,)] = this_result.params[this_conc]
            data['alpha_{}_lb'.format(this_conc,)] = this_result.conf_int(0.05)[0][this_conc]
            data['alpha_{}_ub'.format(this_conc,)] = this_result.conf_int(0.05)[1][this_conc]
            data['alpha_{}_p'.format(this_conc,)] = this_result.pvalues[this_conc]
        else:
            data['alpha_{}'.format(this_conc,)] = .0
            data['alpha_{}_lb'.format(this_conc,)] = .0
            data['alpha_{}_ub'.format(this_conc,)] = .0
            data['alpha_{}_p'.format(this_conc,)] = this_result.pvalues[this_conc]

    return data




