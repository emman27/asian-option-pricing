import math
import numpy
from fixed_call import FixedCall

class AsianVecerFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx * 2, numt, r, sigma, initial_price, strike)
        # We choose maxx to be = q(0)
        self.dx = (2 * self.maxx) / self.numx
        self.set_boundary_conditions()
        self.set_strange_boundary()
        self.j0 = int((self.xi_initial + self.maxx) / self.dx)
        self.xi_min = -self.maxx

    def a(self, t):
        return 0

    def b(self, t):
        return 1

    def initial_value_at_height(self, row):
        return max(-self.maxx + row * self.dx, 0)

    def set_strange_boundary(self):
        for row in range(self.numx + 1):
            for col in range(1, self.numt + 1):
                self.grid.itemset((row, col), self.strange_boundary_value())

    def strange_boundary_value(self):
        return self.q(0)

    def alpha(self, height, time):
        return .5 * self.sigma ** 2 * (
            (1 - math.exp(-self.r * (time + 0.5) * self.dt)) / (self.r * self.maxt)
            + self.maxx - height * self.dx
        ) ** 2 / (2 * self.dx ** 2) * self.dt

    def A_matrix(self, time):
        return super().A_matrix(time, lambda a, b: a, lambda a, b: 1 - 2*a, lambda a, b: a)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a, lambda a, b: -2*a, lambda a, b: a)

    def solve(self):
        return super().solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))

# if __name__ == '__main__':
#     # Constants
#     numx = 200
#     numt = 400
#     maxt = 1
#     r = 0.09
#     s0 = 100

#     sigma = 0.05
#     print('Expected: 13.38, Actual: ' + str(AsianVecerFixedCall(maxt, numx, numt, r, sigma, s0, 90).solve()))
#     print('Expected: 8.81, Actual: ' + str(AsianVecerFixedCall(maxt, numx, numt, r, sigma, s0, 95).solve()))
#     print('Expected: 4.33, Actual: ' + str(AsianVecerFixedCall(maxt, numx, numt, r, sigma, s0, 100).solve()))
#     print('Expected: 0.88, Actual: ' + str(AsianVecerFixedCall(maxt, numx, numt, r, sigma, s0, 105).solve()))
#     print('Expected: 0.06, Actual: ' + str(AsianVecerFixedCall(maxt, numx, numt, r, sigma, s0, 110).solve()))
