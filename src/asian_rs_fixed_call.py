from fixed_call import FixedCall
import math
import numpy

class AsianRSFixedCall(FixedCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.dx = self.maxx / self.numx
        self.j0 = int(self.xi_initial / self.dx)
        self.set_boundary_conditions()

    def initial_value_at_bottom(self, time):
        top = (1 - math.exp(-self.r * time * self.dt))
        bottom = (self.r * self.maxt)
        return top / bottom

    def a(self, t):
        return -self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    # Methods to calculate the abstracted alpha and beta variables. See report.
    def alpha(self, height, time):
        return .25 * self.sigma**2 * height**2 * self.dt

    def beta(self, height, time):
        return (height * self.r * self.dx + 1 / self.maxt) * self.dt  / (4 * self.dx)

    # The relation between two columns are given by
    # (I-B) * R = A * L
    # where they follow the original relation
    # R = AL + BR. See report.
    # Here time is set to 0 for the super() call since time is irrelevant for this option
    def A_matrix(self):
        return super().A_matrix(0, lambda a, b: a + b, lambda a, b: 1 - 2*a, lambda a, b: a - b)

    def B_matrix(self):
        return super().B_matrix(0, lambda a, b: a + b, lambda a, b: -2 * a, lambda a, b: a - b)

    def solve(self):
        a_mat = self.A_matrix()
        b_mat = self.B_matrix()
        self.fake_super_solve(lambda time: numpy.identity(self.numx + 1) - b_mat, lambda time: a_mat)
        return self.interpolate()

    def fake_super_solve(self, left_multiplier, right_multiplier):
        for col in range(self.numt):
            l = self.grid[:, col]
            r = self.grid[:, col + 1]
            new = numpy.linalg.solve(left_multiplier(col), right_multiplier(col) * l)
            for row in range(1, self.numx + 1):
                self.grid.itemset((row, col + 1), new[row])
