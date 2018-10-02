import numpy as np
import statsmodels
from pandas import  DataFrame

def geek_regression(df,
                    concentrations,
                    reference_concentrations,
                    rate_constant,
                    verbose=True):
    """

    :param df: dataframe with measured rate constant data
    :param concentrations: list of column names of the concentrations
    :param reference_concentrations: numpy array of reference concentrations
    :param rate_constant: column name of the rate constant data
    :param verbose: True or False print regression output
    :return: A dict with the geek parameters
    """

    if len(concentrations) != len(reference_concentrations):
        raise ValueError("You must provide a reference concentration for each concentration column")

    # Prepare multivariate linear regression
    X = np.log(df[concentrations]/reference_concentrations)
    Y = df[rate_constant]
    X = statsmodels.add_constant(X)

    # Get the odinary lineast squares residuals
    ols_result = statsmodels.OLS(Y, X).fit()
    ols_resid = DataFrame(dict(fittedvalues=ols_result.fittedvalues,
                               residuals=ols_result.resid))

    # Group the values by their residual with in relative accuacy of 0.1
    ols_resid['value_groups'] = (ols_resid['fittedvalues'] / max(abs(ols_resid['fittedvalues']))).round(1)

    # Calculate the conditional standard deviation of the fitted values
    conditional_std = ols_resid.groupby('value_groups')['residuals'].std()
    ols_resid['std'] = np.repeat(np.NaN, len(ols_resid))
    for this_group in conditional_std.index.values:
        ols_resid['std'][ols_resid['value_groups'] == this_group] = conditional_std.ix[this_group]

    # Fit with weights antiproportinal to the standart deviation
    weights = 1.0 / ols_resid['std']

    # Fit the model
    this_lin_model = statsmodels.WLS(Y, X, weights=weights)
    this_result = this_lin_model.fit()

    if verbose:
        print(this_result.summary())

    # Data for the
    data = {'rate_constant': rate_constant,
            'beta': this_result.params['const'],
            'lb_beta': this_result.conf_int(0.05)[0]['const'],
            'ub_beta': this_result.conf_int(0.05)[1]['const'],
            'p_beta': this_result.pvalues['const'],
            }

    for this_conc in concentrations:
        data['alpha_{}_{}'.format(this_conc,rate_constant)] = this_result.params[this_conc]
        data['lb_alpha_{}_{}'.format(this_conc, rate_constant)] = this_result.conf_int(0.05)[0][this_conc]
        data['ub_alpha_{}_{}'.format(this_conc, rate_constant)] = this_result.conf_int(0.05)[1][this_conc]
        data['p_alpha_{}_{}'.format(this_conc, rate_constant)] = this_result.pvalues[this_conc]

    return data




