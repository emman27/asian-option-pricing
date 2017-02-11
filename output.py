#!/usr/bin/python3
import csv
import progressbar
from asian_rs_fixed_call import AsianRSFixedCall
from asian_vecer_fixed_call import AsianVecerFixedCall
from asian_dl_fixed_call import AsianDLFixedCall
from asian_new_fixed_call import AsianNewFixedCall

from asian_rs_float_call import AsianRSFloatCall
from asian_vecer_float_call import AsianVecerFloatCall
from asian_dl_float_call import AsianDLFloatCall
from asian_new_float_call import AsianNewFloatCall

if __name__ == "__main__":
    # Fixed Strike Constants
    NUMX = 200
    NUMT = 400
    MAXT = 1
    R = 0.09
    S0 = 100
    SIGMAS = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    STRIKES = [90, 95, 100, 105, 110]
    METHODS = [AsianRSFixedCall, AsianVecerFixedCall, AsianDLFixedCall, AsianNewFixedCall]
    BENCHMARKS = [
        [13.38, 8.81, 4.31, 0.96, 0.05],
        [13.39, 8.91, 4.92, 2.07, 0.63],
        [13.83, 10.00, 6.78, 4.30, 2.55],
        [14.98, 11.66, 8.83, 6.52, 4.70],
        [16.50, 13.51, 10.92, 8.73, 6.90],
        [18.19, 15.44, 13.03, 10.93, 9.12],
        [19.96, 17.41, 15.13, 13.11, 11.34]
    ]
    HEADERS = ['Volatility', 'Strike', 'Benchmark', 'Rogers-Shi', 'Vecer', 'Dubois Lelievre', 'New']


    # FLoating strike constants
    R_FLOAT = 0.1
    S0_FLOAT = 100
    SIGMAS_FLOAT = [0.3, 0.5]
    T0S = [0.1, 0.3, 0.5, 0.7, 0.9]
    CURRENT_AVERAGES = [90, 100, 110]
    METHODS_FLOAT = [AsianRSFloatCall, AsianVecerFloatCall, AsianDLFloatCall, AsianNewFloatCall]
    HEADERS_FLOAT = ['Volatility', 't0', 'Average', 'Benchmark', 'Rogers-Shi', 'Vecer', 'Dubois Lelievre', 'New']
    BENCHMARKS_FLOAT = [
        [
            [9.86, 9.35, 8.85],
            [10.69, 9.06, 7.61],
            [11.20, 8.31, 6.00],
            [11.19, 6.86, 3.86],
            [10.38, 4.07, 1.04]
        ],# sigma = 0.3
        [
            [14.09, 13.65, 13.23],
            [14.71, 13.33, 12.08],
            [14.88, 12.42, 10.32],
            [14.21, 10.49, 7.58],
            [11.91, 6.44, 3.05]
        ] # sigma = 0.5
    ]

    output = []
    output.append(HEADERS)
    # len(SIGMAS)*len(STRIKES)*len(METHODS) +
    bar = progressbar.ProgressBar(max_value = len(SIGMAS_FLOAT)*len(T0S)*len(CURRENT_AVERAGES)*len(METHODS_FLOAT))
    k = 0
    for j in range(len(SIGMAS)):
        sigma = SIGMAS[j]
        for i in range(len(STRIKES)):
            strike = STRIKES[i]
            results = [sigma, strike]
            results.append(BENCHMARKS[j][i])
            for method in METHODS:
                results.append(round(method(MAXT, NUMX, NUMT, R, sigma, S0, strike).solve(), 2))
                k += 1
                bar.update(k)
            output.append(results)

    with open('output_fixed.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        for row in output:
            writer.writerow(row)

    # Floating Strike output
    output = []
    output.append(HEADERS_FLOAT)
    for j in range(len(SIGMAS_FLOAT)):
        sigma = SIGMAS_FLOAT[j]
        for n in range(len(T0S)):
            t0 = T0S[n]
            for i in range(len(CURRENT_AVERAGES)):
                avr = CURRENT_AVERAGES[i]
                results = [sigma, t0, avr]
                results.append(BENCHMARKS_FLOAT[j][n][i])
                for method in METHODS_FLOAT:
                    results.append(round(method(MAXT - t0, NUMX, NUMT, R_FLOAT, sigma, S0, avr, t0).solve(), 2))
                    k += 1
                    bar.update(k)
                output.append(results)

    with open('output_floating.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        for row in output:
            writer.writerow(row)

