""" Obtain results from csv and extrapolate to CBS & FCI limit."""
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm

np.set_printoptions(precision=12)


def printCorrelationEnergy(statsResult):
    coefs = statsResult.params.values
    stdevs = statsResult.bse
    print('Correlation Energy: ' + str(coefs[0]) + ' +- ' + str(stdevs[0]))


def BEWRegression(X, y, e, title):
    augX = sm.add_constant(X)

    # Backward elimination.
    print('\n' + '#' * 80)
    print('REG: ' + title)
    print('-' * 80)
    print('Backward elimination:')
    results = sm.WLS(y, augX, weights=1.0 / np.square(e)).fit()
    intercept = results.params.values[0]
    variance = np.square(
        np.dot(augX, np.abs(results.params.values)) + intercept) + np.square(e)
    iteration = 0
    while True:
        results = sm.WLS(y, augX, weights=1.0 / variance).fit()
        printCorrelationEnergy(results)
        intercept = results.params.values[0]
        variance = np.square(
            np.dot(augX, np.abs(results.params.values)) + intercept) + np.square(e)

        iteration = iteration + 1
        if iteration < 5:
            continue

        maxPIndex = np.argmax(results.pvalues)
        maxP = results.pvalues[maxPIndex]

        if maxP < 0.5:
            break
        print('Eliminate: ' + maxPIndex)
        print('P > |t|: ' + str(maxP))
        augX.drop(maxPIndex, axis=1, inplace=True)

    # Weighted OLS
    print('\n[FINAL Weighted OLS]')
    variance = np.square(
        np.dot(augX, np.abs(results.params.values)) + intercept) + np.square(e)
    results = sm.WLS(y, augX, weights=1.0 / variance).fit()
    # print(results.summary())
    printCorrelationEnergy(results)


def main():
    """main function"""
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

    # Add cross terms.
    selectedParameters = parameters[:]
    for i in range(len(parameters)):
        for j in range(i, len(parameters)):
            for k in range(j, len(parameters)):
                column = parameters[i] + ' * ' + \
                    parameters[j] + ' * ' + parameters[k]
                selectedParameters.append(column)
                data[column] = data[parameters[i]] * \
                    data[parameters[j]] * data[parameters[k]]

    # Estimate intercept.
    X = data[selectedParameters]
    y = data['energy_corr']
    e = data['uncert']

    BEWRegression(X, y, e, 'all data')

    for i, parameter in enumerate(parameters):
        maxValue = X.max()[i]
        keep = X[parameter] != maxValue
        X_rmax = X[keep]
        y_rmax = y[keep]
        e_rmax = e[keep]
    BEWRegression(X_rmax, y_rmax, e_rmax, 'accu data')

    for i, parameter in enumerate(parameters):
        minValue = X.min()[i]
        keep = X[parameter] != minValue
        X_rmin = X[keep]
        y_rmin = y[keep]
        e_rmin = e[keep]
    BEWRegression(X_rmin, y_rmin, e_rmin, 'verify data')

    for i, parameter in enumerate(parameters):
        maxValue = X_rmin.max()[i]
        keep = X_rmin[parameter] != maxValue
        X_rminmax = X_rmin[keep]
        y_rminmax = y_rmin[keep]
        e_rminmax = e_rmin[keep]
    BEWRegression(X_rminmax, y_rminmax, e_rminmax, 'verify accu data')


if __name__ == '__main__':
    main()
