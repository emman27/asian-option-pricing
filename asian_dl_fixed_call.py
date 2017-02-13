from fixed_call import FixedCall
import math
import numpy

class AsianDLFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.xi_initial = self.xi(self.s0, 0)
        self.j0 = round(numx/3)
        self.dx = self.xi_initial / self.j0
        self.set_boundary_conditions()

    def a(self, t):
        return t/self.maxt - self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    def set_initial_boundary(self):
        for row in range(self.numx):
            self.grid.itemset((row, 0), self.initial_value_at_height(row))

    def initial_value_at_height(self, row):
        return max(1 - row * self.dx, 0)

    def set_bottom_boundary(self):
        for col in range(self.numt + 1):
            self.grid.itemset((0, col), self.initial_value_at_time(col))

    def initial_value_at_time(self, col):
        return (
            (1 - math.exp(-self.r * col * self.dt)) / (self.r * self.maxt) +
            math.exp(-self.r * col * self.dt) * (self.maxt - col * self.dt) / self.maxt
        )

    def alpha(self, row, col):
        return 0.5 * self.sigma ** 2 * (1 - row * self.dx - ((col + .5) * self.dt) / self.maxt) ** 2 / (2 * self.dx ** 2) * self.dt

    def beta(self, row, col):
        return self.r * self.dt * (1 - row * self.dx - ((col + .5) * self.dt) / self.maxt) / (4 * self.dx)

    def A_matrix(self, time):
        A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i, time)
            b = self.beta(i, time)
            if i - 1 >= 0:
                A.itemset((i, i-1), a - b)
            A.itemset((i, i), 1 - 2 * a)
            if i + 1 < self.numx:
                A.itemset((i, i + 1), a + b)
        return A

    def B_matrix(self, time):
        B = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i, time)
            b = self.beta(i, time)
            if i - 1 >= 0:
                B.itemset((i, i-1), a - b)
            B.itemset((i, i), -2 * a)
            if i + 1 < self.numx:
                B.itemset((i, i + 1), a + b)
        return B

    def solve(self):
        super().solve(lambda time: numpy.identity(self.numx) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.s0 * self.grid[self.j0, self.numt]

if __name__ == '__main__':
    # Constants
    numx = 200
    numt = 400
    maxt = 1
    r = 0.09
    s0 = 100

    sigma = 0.05
    # print('Expected: 13.38, Actual: ' + str(AsianDLFixedCall(maxt, numx, numt, r, sigma, s0, 90).solve()))
    # print('Expected: 8.81, Actual: ' + str(AsianDLFixedCall(maxt, numx, numt, r, sigma, s0, 95).solve()))
    # print('Expected: 4.22, Actual: ' + str(AsianDLFixedCall(maxt, numx, numt, r, sigma, s0, 100).solve()))
    # print('Expected: 1.00, Actual: ' + str(AsianDLFixedCall(maxt, numx, numt, r, sigma, s0, 105).solve()))
    # print('Expected: 0.09, Actual: ' + str(AsianDLFixedCall(maxt, numx, numt, r, sigma, s0, 110).solve()))
