#!/usr/bin/python3
import numpy
import math

# Constants
numx = 200
numt = 400
maxt = 1
r = 0.09
s0 = 100

def solve(numx, numt, maxt, r, sigma, strike, s0):
    # Helpers
    def initial_zero_value_at_time(time):
        '''
        Calculates the zero-height value of the FDS grid at given time @time
        '''
        top = (1 - math.exp(-r * time))
        bottom = (r * maxt)
        return top / bottom

    def alpha(height):
        return .25 * sigma**2 * height**2 * dt

    def beta(height):
        return (height * r * dx + 1 / maxt) * dt  / (4 * dx)

    def A_matrix():
        A = numpy.matrix([[0] * numx] * numx, dtype = numpy.float64)
        for i in range(numx):
            a = alpha(i)
            b = beta(i)
            try:
                A.itemset((i, i - 1), (a + b))
                A.itemset((i, i), (1 - 2*a))
                A.itemset((i, i + 1), (a - b))
            except Exception:
                pass
        return A

    def B_matrix():
        B = numpy.matrix([[0] * numx] * numx, dtype = numpy.float64)
        for i in range(numx):
            a = alpha(i)
            b = beta(i)
            try:
                B.itemset((i, i - 1), (a + b))
                B.itemset((i, i), (-2 * a))
                B.itemset((i, i + 1), (a - b))
            except:
                pass
        return B

    def xi(s, t):
        return a(t) + b(t) * (avr(t) - strike * math.exp(-r * (maxt - t))) / s

    def avr(t):
        return q(t) * s0

    def a(t):
        return -q(t) * b(t)

    def b(t):
        return -math.exp(r*(maxt - t))

    def q(t):
        return (1 - math.exp(-r*(maxt - t))) / (r * maxt)

    dt = maxt / numt
    xi_initial = xi(s0, 0)
    # Empirical testing shows that a j value of 3 makes a good choice.
    j0 = numx // 3
    dx = strike / s0 / j0
    a_mat = A_matrix()
    b_mat = B_matrix()

    # Setup
    # There are in fact n+1 time points to account for the fact that there needs to be a state for t = 0
    grid = numpy.matrix([[0] * (numt + 1)] * numx, dtype=numpy.float64)

    for col in range(numt+1):
        grid.itemset((0, col), initial_zero_value_at_time((col) * dt))

    # Solving
    for col in range(numt):
        l = grid[:, col]
        r = grid[:, col + 1]
        new = numpy.linalg.solve(numpy.identity(numx) - b_mat, a_mat * l)
        for row in range(1, numx):
            grid.itemset((row, col + 1), new[row])

    return (s0 * grid[j0, numt])

print('With σ = 0.05')
print('Expected: 13.07, Actual: ' + str(solve(numx, numt, maxt, r, 0.05, 90, s0)))
print('Expected: 7.82, Actual: ' + str(solve(numx, numt, maxt, r, 0.05, 95, s0)))
print('Expected: 3.91, Actual: ' + str(solve(numx, numt, maxt, r, 0.05, 100, s0)))
print('Expected: 1.65, Actual: ' + str(solve(numx, numt, maxt, r, 0.05, 105, s0)))
print('Expected: 0.59, Actual: ' + str(solve(numx, numt, maxt, r, 0.05, 110, s0)))

print('With σ = 0.1')
print('Expected: 13.11, Actual: ' + str(solve(numx, numt, maxt, r, 0.1, 90, s0)))
print('Expected: 8.37, Actual: ' + str(solve(numx, numt, maxt, r, 0.1, 95, s0)))
print('Expected: 4.68, Actual: ' + str(solve(numx, numt, maxt, r, 0.1, 100, s0)))
print('Expected: 2.30, Actual: ' + str(solve(numx, numt, maxt, r, 0.1, 105, s0)))
print('Expected: 1, Actual: ' + str(solve(numx, numt, maxt, r, 0.1, 110, s0)))
