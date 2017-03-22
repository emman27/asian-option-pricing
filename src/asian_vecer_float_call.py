from floating_call import FloatingCall
import numpy
import math

class AsianVecerFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx * 2, numt, r, sigma, initial_price, old_average, t0)
        self.dx = self.maxx * 2 / self.numx
        self.j0 = int((self.xi_initial + self.maxx) / self.dx)
        self.xi_min = -self.maxx
        self.set_boundary_conditions()

    def a(self, t):
        return 0

    def b(self, t):
        return 1

    def initial_value_at_height(self, row):
        return max(1 + self.maxx - row * self.dx, 0)

    def initial_value_at_bottom(self, col):
        return 1 + self.maxx

    def alpha(self, height, time):
        return 1/2 * self.sigma**2 * ((1 - math.exp(-self.r*(time + 0.5) * self.dt))/(self.r * (self.t0 + self.maxt)) + (self.numx / 2 - height) * self.dx)**2 / (2 * self.dx**2) * self.dt

    def A_matrix(self, time):
        return super().A_matrix(time, lambda a, b: a, lambda a, b: 1 - 2 * a, lambda a, b: a)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a, lambda a, b: -2 * a, lambda a, b: a)

    def solve(self):
        return super().solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))
