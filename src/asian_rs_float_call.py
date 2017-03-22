from floating_call import FloatingCall
import numpy
import math

class AsianRSFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx * 2, numt, r, sigma, initial_price, old_average, t0)
        self.dx = self.maxx * 2 / self.numx
        self.j0 = int((self.xi_initial + self.maxx) / self.dx)
        self.xi_min = -self.maxx
        self.set_boundary_conditions()

    def a(self, t):
        return -self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    def initial_value_at_top(self, col):
        return 1 - (1 - math.exp(-self.r * col * self.dt)) / (self.r * (self.t0 + self.maxt)) + math.exp(-self.r * col * self.dt) * self.maxx

    def initial_value_at_height(self, row):
        return max(1 - self.maxx + row * self.dx, 0)

    def alpha(self, height, time):
        return 1/2 * self.sigma**2 * (-self.numx / 2 + height)**2 * self.dx**2 / (2*self.dx**2) * self.dt

    def beta(self, height, time):
        return -( self.r * (-self.numx / 2 + height) * self.dx + 1 / (self.t0 + self.maxt) ) / (4 * self.dx) * self.dt

    def A_matrix(self):
        return super().A_matrix(0, lambda a, b: a - b, lambda a, b: 1 - 2*a, lambda a, b: a + b)

    def B_matrix(self):
        return super().B_matrix(0, lambda a, b: a - b, lambda a, b: -2 * a, lambda a, b: a + b)

    def solve(self):
        a_mat = self.A_matrix()
        b_mat = numpy.identity(self.numx + 1) - self.B_matrix()
        return super().solve(lambda time: b_mat, lambda time: a_mat)
