""" Obtain results from csv and extrapolate to CBS & FCI limit."""
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
import itertools

np.set_printoptions(precision=12)


def extrapolate(x, y, e):
    augX = x.copy()
    augX = sm.add_constant(augX)
    augX['x^2'] = np.power(augX[augX.columns[1]], 2.0)
    # Initial fit.
    results = sm.WLS(y, augX, weights=1.0 / np.square(e)).fit()
    intercept = results.params.values[0]
    variance = np.square(
        np.dot(augX, np.abs(results.params.values)) + intercept) + np.square(e)

    for i in range(5):
        results = sm.WLS(y, augX, weights=1.0 / variance).fit()
        intercept = results.params.values[0]
        variance = np.square(
            np.dot(augX, np.abs(results.params.values)) + intercept) + np.square(e)

    intercept = results.params.values[0]
    uncert = results.bse[0]
    return [intercept, uncert]


def extrapolate_to(parameter, ref_parameters, data):
    print('#' * 80)
    print('Extrapolating: ' + parameter)
    print('Ref: ' + str(ref_parameters))
    print()

    unique_values_all = []
    for ref_parameter in ref_parameters:
        print('Parameter: ' + ref_parameter)
        unique_values = data[ref_parameter].unique()
        print(unique_values)
        unique_values_all.append(unique_values.tolist())

    res_columns = ref_parameters + ['energy_corr', 'uncert']
    new_data = pd.DataFrame(
        data=None, columns=res_columns, dtype=np.float64)

    print(unique_values_all)
    categories = list(itertools.product(*unique_values_all))

    for category in categories:
        locs = (data[parameter] != np.NaN)
        for index, ref_parameter in enumerate(ref_parameters):
            locs = locs & (data[ref_parameter] == category[index])
        category_data = data[locs]
        category_extrapolate = extrapolate(
            category_data[parameter], category_data['energy_corr'], category_data['uncert'])
        new_data.loc[len(new_data)] = list(category) + category_extrapolate

    if len(categories) == 0:
        category_data = data
        category_extrapolate = extrapolate(
            category_data[parameter], category_data['energy_corr'], category_data['uncert'])
        new_data.loc[len(new_data)] = category_extrapolate
    print(new_data)
    return new_data


def main():
    # Check and read res file.
    res_file = 'pt_result.csv'
    if len(sys.argv) == 2:
        res_file = sys.argv[1]

    parameters = ['n_orbs_var_inv', 'eps_var', 'n_orbs_pt_inv', 'eps_pt']

    # Read raw data.
    data = pd.read_csv(res_file)

    # Add inverse terms.
    data['n_orbs_var_inv'] = 1.0 / data['n_orbs_var']
    data['n_orbs_pt_inv'] = 1.0 / data['n_orbs_pt']

    # Remove parameters not enough for extrapolation.
    for parameter in parameters:
        if data[parameter].value_counts().size < 3:
            parameters.remove(parameter)

    ref_parameters = parameters.copy()
    data_extrapolate = data.copy()
    for parameter in parameters[::-1]:
        ref_parameters.remove(parameter)
        data_extrapolate = extrapolate_to(
            parameter, ref_parameters, data_extrapolate)

    print('#' * 80)
    print('#' * 80)
    print('#' * 80)

    for i, parameter in enumerate(parameters):
        maxValue = data[parameter].max()
        keep = data[parameter] != maxValue
        data = data[keep]
    ref_parameters = parameters.copy()
    data_extrapolate = data.copy()
    for parameter in parameters[::-1]:
        ref_parameters.remove(parameter)
        data_extrapolate = extrapolate_to(
            parameter, ref_parameters, data_extrapolate)


if __name__ == '__main__':
    main()
