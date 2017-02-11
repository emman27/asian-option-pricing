from fixed_call import FixedCall
import math
import numpy

class AsianRSFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.j0 = round(numx/3)
        self.xi_initial = self.xi(self.s0, 0)
        self.dx = self.strike / self.s0 / self.j0
        self.set_boundary_conditions()

    def set_bottom_boundary(self):
        '''
        Sets the boundary values for the FDS grid when the height is 0
        '''
        for col in range(self.numt + 1):
            time = col * self.dt
            self.grid.itemset((0, col), self.initial_zero_value_at_time(time))

    def initial_zero_value_at_time(self, time):
        '''
        Calculates the zero-height value of the FDS grid at given time @time
        '''
        top = (1 - math.exp(-self.r * time))
        bottom = (self.r * self.maxt)
        return top / bottom

    def a(self, t):
        return -self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    # Methods to calculate the abstracted alpha and beta variables. See report.
    def alpha(self, height):
        return .25 * self.sigma**2 * height**2 * self.dt

    def beta(self, height):
        return (height * self.r * self.dx + 1 / self.maxt) * self.dt  / (4 * self.dx)

    # The relation between two columns are given by
    # (I-B) * R = A * L
    # where they follow the original relation
    # R = AL + BR. See report.
    def A_matrix(self):
        A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i)
            b = self.beta(i)
            if i - 1 >= 0:
                A.itemset((i, i - 1), (a + b))
            A.itemset((i, i), (1 - 2*a))
            if i + 1 < self.numx:
                A.itemset((i, i + 1), (a - b))
        return A

    def B_matrix(self):
        B = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i)
            b = self.beta(i)
            if i - 1 >= 0:
                B.itemset((i, i - 1), (a + b))
            B.itemset((i, i), (-2 * a))
            if i + 1 < self.numx:
                B.itemset((i, i + 1), (a - b))
        return B

    def solve(self):
        a_mat = self.A_matrix()
        b_mat = self.B_matrix()
        self.fake_super_solve(lambda time: numpy.identity(self.numx) - b_mat, lambda time: a_mat)
        return self.s0 * self.grid[self.j0, self.numt]

    def fake_super_solve(self, left_multiplier, right_multiplier):
        for col in range(self.numt):
            l = self.grid[:, col]
            r = self.grid[:, col + 1]
            new = numpy.linalg.solve(left_multiplier(col), right_multiplier(col) * l)
            for row in range(1, self.numx):
                self.grid.itemset((row, col + 1), new[row])
# Constants
numx = 200
numt = 400
maxt = 1
r = 0.09
s0 = 100

sigma = 0.05
# print('Expected: 13.07, Actual: ' + str(AsianRSFixedCall(maxt, numx, numt, r, sigma, s0, 90).solve()))
# print('Expected: 7.82, Actual: ' + str(AsianRSFixedCall(maxt, numx, numt, r, sigma, s0, 95).solve()))
# print('Expected: 3.91, Actual: ' + str(AsianRSFixedCall(maxt, numx, numt, r, sigma, s0, 100).solve()))
# print('Expected: 1.65, Actual: ' + str(AsianRSFixedCall(maxt, numx, numt, r, sigma, s0, 105).solve()))
# print('Expected: 0.59, Actual: ' + str(AsianRSFixedCall(maxt, numx, numt, r, sigma, s0, 110).solve()))
