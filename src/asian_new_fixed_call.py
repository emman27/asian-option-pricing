from fixed_call import FixedCall
import math
import numpy

class AsianNewFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.maxx = abs(5 * self.xi_initial)
        self.dx = self.maxx / self.numx
        self.j0 = self.numx + round(self.xi_initial / self.dx)
        self.set_boundary_conditions()

    def a(self, t):
        return -self.q(t)

    def b(self, t):
        return 1

    def initial_value_at_top(self, col):
        return (1 - math.exp(-self.r * col * self.dt)) / (self.r * self.maxt)

    def alpha(self, row, col):
        return .5 * self.sigma**2 * (-self.maxx + row * self.dx) ** 2 / (2 * self.dx**2) * self.dt

    def beta(self, row, col):
        return math.exp(-self.r * (col + 0.5) * self.dt) / (self.maxt * 4 * self.dx) * self.dt

    def A_matrix(self, time):
        return super().A_matrix(time, lambda a, b: a - b, lambda a, b: 1 - 2*a, lambda a, b: a + b)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a - b, lambda a, b: - 2*a, lambda a, b: a + b)

    def solve(self):
        self.fake_super_solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.s0 * self.grid[self.j0, self.numt]

    def fake_super_solve(self, left_multiplier, right_multiplier):
        for col in range(self.numt):
            l = self.grid[:, col]
            r = self.grid[:, col + 1]
            new = numpy.linalg.solve(left_multiplier(col), right_multiplier(col) * l)
            for row in range(self.numx):
                self.grid.itemset((row, col + 1), new[row])

# if __name__ == "__main__":
#     # Constants
#     numx = 200
#     numt = 400
#     maxt = 1
#     r = 0.09
#     s0 = 100

#     sigma = 0.05
#     print('Expected: 13.07, Actual: ' + str(AsianNewFixedCall(maxt, numx, numt, r, sigma, s0, 90).solve()))
#     print('Expected: 7.81, Actual: ' + str(AsianNewFixedCall(maxt, numx, numt, r, sigma, s0, 95).solve()))
#     print('Expected: 3.91, Actual: ' + str(AsianNewFixedCall(maxt, numx, numt, r, sigma, s0, 100).solve()))
#     print('Expected: 1.65, Actual: ' + str(AsianNewFixedCall(maxt, numx, numt, r, sigma, s0, 105).solve()))
#     print('Expected: 0.59, Actual: ' + str(AsianNewFixedCall(maxt, numx, numt, r, sigma, s0, 110).solve()))
