from fixed_call import FixedCall
import math
import numpy

class AsianNewVecerFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.xi_initial = self.xi(self.s0, 0)
        self.j0 = round(self.numx/3)
        self.dx = self.xi_initial / self.j0
        self.maxx = self.dx * (self.numx - 1)
        self.set_boundary_conditions()

    def xi(self, s, t):
        return math.exp(-self.r * (self.maxt - t)) * s / self.avr(s, t)

    def q(self, t):
        return (math.exp(self.r * self.maxt - t) - 1) / (self.r * self.maxt)

    def avr(self, s, t):
        return self.q(t) * math.exp(-self.r * (self.maxt - t)) * s

    def set_initial_boundary(self):
        for row in range(self.numx):
            self.grid.itemset((row, 0), self.initial_value_at_height(row))

    def initial_value_at_height(self, row):
        return max(1 - self.strike * row * self.dx / self.s0, 0)

    def set_bottom_boundary(self):
        for col in range(self.numt):
            self.grid.itemset((0, col), self.initial_value_at_time(col))

    def initial_value_at_time(self, col):
        return 1

    def alpha(self, height, time):
        return .5 * self.sigma ** 2 * (height * self.dx) ** 2 * (
            (math.exp(self.r * (time + 0.5) * self.dt) - 1) /
            (self.r * self.maxt) *
            (height * self.dx) -
            1
        )**2 / (2 * self.dx**2) * self.dt

    def A_matrix(self, time):
        A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i, time)
            if i - 1 >= 0:
                A.itemset((i, i-1), a)
            A.itemset((i, i), 1 - 2 * a)
            if i + 1 < self.numx:
                A.itemset((i, i + 1), a)
        return A

    def B_matrix(self, time):
        B = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i, time)
            if i - 1 >= 0:
                B.itemset((i, i-1), a)
            B.itemset((i, i), -2 * a)
            if i + 1 < self.numx:
                B.itemset((i, i + 1), a)
        return B

    def solve(self):
        super().solve(lambda time: numpy.identity(self.numx) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.avr(self.s0, 0) * self.grid[self.j0, self.numt]

if __name__ == '__main__':
    # Constants
    numx = 200
    numt = 400
    maxt = 1
    r = 0.09
    s0 = 100

    sigma = 0.05
    # print('Expected: 13.07, Actual: ' + str(AsianNewVecerFixedCall(maxt, numx, numt, r, sigma, s0, 90).solve()))
    # print('Expected: 7.81, Actual: ' + str(AsianNewVecerFixedCall(maxt, numx, numt, r, sigma, s0, 95).solve()))
    # print('Expected: 3.91, Actual: ' + str(AsianNewVecerFixedCall(maxt, numx, numt, r, sigma, s0, 100).solve()))
    # print('Expected: 1.65, Actual: ' + str(AsianNewVecerFixedCall(maxt, numx, numt, r, sigma, s0, 105).solve()))
    # print('Expected: 0.59, Actual: ' + str(AsianNewVecerFixedCall(maxt, numx, numt, r, sigma, s0, 110).solve()))
