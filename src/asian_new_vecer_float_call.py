from floating_call import FloatingCall
import math
import numpy

class AsianNewVecerFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, old_average, t0)
        self.dx = self.maxx / self.numx
        self.j0 = int(self.xi_initial / self.dx)
        self.set_boundary_conditions()

    def xi(self, s, t):
        return s / self.avr(t)

    def initial_value_at_height(self, row):
        return max(row * self.dx - 1, 0)

    def initial_value_at_top(self, col):
        return self.maxx - 1

    def alpha(self, height, time):
        return .5 * self.sigma**2 * (height * self.dx)**2 * (
            (1 - math.exp(-self.r * (time + 0.5) * self.dt)) /
            (self.r * (self.t0 + self.maxt)) *
            height * self.dx -
            1
        ) ** 2 * self.dt / (2 * self.dx**2)

    def A_matrix(self, time):
        return super().A_matrix(time, lambda a, b: a, lambda a, b: 1 - 2 * a, lambda a, b: a)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a, lambda a, b: - 2 * a, lambda a, b: a)

    def solve(self):
        self.fake_super_solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.interpolate()

    def fake_super_solve(self, left_multiplier, right_multiplier):
        for col in range(self.numt):
            l = self.grid[:, col]
            r = self.grid[:, col + 1]
            new = numpy.linalg.solve(left_multiplier(col), right_multiplier(col) * l)
            for row in range(self.numx):
                self.grid.itemset((row, col + 1), new[row])
