from fixed_call import FixedCall
import math
import numpy

class AsianDLFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.dx = self.maxx / self.numx
        self.j0 = round(self.xi_initial / self.dx)
        self.set_boundary_conditions()

    def a(self, t):
        return t/self.maxt - self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    def initial_value_at_height(self, row):
        return max(1 - row * self.dx, 0)

    def initial_value_at_bottom(self, col):
        return (
            (1 - math.exp(-self.r * col * self.dt)) / (self.r * self.maxt) +
            math.exp(-self.r * col * self.dt) * (self.maxt - col * self.dt) / self.maxt
        )

    def alpha(self, row, col):
        return 0.5 * self.sigma ** 2 * (1 - row * self.dx - ((col + .5) * self.dt) / self.maxt) ** 2 / (2 * self.dx ** 2) * self.dt

    def beta(self, row, col):
        return self.r * self.dt * (1 - row * self.dx - ((col + .5) * self.dt) / self.maxt) / (4 * self.dx)

    def A_matrix(self, time):
        return super().A_matrix(time, lambda a, b: a - b, lambda a, b: 1 - 2*a, lambda a, b: a + b)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a - b, lambda a, b: - 2*a, lambda a, b: a + b)

    def solve(self):
        super().solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.s0 * self.grid[self.j0, self.numt]
