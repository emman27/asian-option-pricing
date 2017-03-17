import numpy
import math

class Option:
    def __init__(self, maxt, numx, numt, r, sigma, s0, avr, t0):
        self.maxt = maxt
        self.numx = numx
        self.numt = numt
        self.r = r
        self.sigma = sigma
        self.s0 = s0
        self.old_average = avr
        self.t0 = t0
        self.dt = maxt / float(numt)
        self.xi_min = 0
        self.xi_initial = self.xi(self.s0, 0)
        self.maxx = 2
        # One extra time point to account for t = 0
        self.grid = numpy.matrix([[0] * (self.numt + 1)] * (self.numx + 1), dtype = numpy.float64)

    def q(self, t):
        return (1 - math.exp(-self.r*(self.maxt - t))) / (self.r * (self.t0 + self.maxt))

    def avr(self, t):
        return self.q(t) * self.s0 + math.exp(-self.r * (self.maxt - t)) * self.t0 * self.old_average / (self.t0 + self.maxt)

    def set_boundary_conditions(self):
        self.set_bottom_boundary()
        self.set_top_boundary()
        self.set_initial_boundary()

    def solve(self, left_multiplier, right_multiplier):
        for col in range(self.numt):
            l = self.grid[:, col]
            r = self.grid[:, col + 1]
            new = numpy.linalg.solve(left_multiplier(col), right_multiplier(col) * l)
            for row in range(self.numx + 1):
                self.grid.itemset((row, col + 1), new[row])
        return self.interpolate()

    def interpolate(self):
        y0 = self.s0 * self.grid[self.j0, self.numt]
        y1 = self.s0 * self.grid[self.j0 + 1, self.numt]
        x = self.xi_initial
        x0 = self.j0 * self.dx + self.xi_min
        x1 = (self.j0 + 1) * self.dx + self.xi_min
        return y0 + (x - x0) * (y1 - y0) / (x1 - x0)

    def set_bottom_boundary(self):
        for i in range(self.numt + 1):
            self.grid.itemset((0, i), self.initial_value_at_bottom(i))

    def set_top_boundary(self):
        for i in range(self.numt + 1):
            self.grid.itemset((self.numx, i), self.initial_value_at_top(i))

    def set_initial_boundary(self):
        for j in range(self.numx + 1):
            self.grid.itemset((j, 0), self.initial_value_at_height(j))

    def A_matrix(self, time, lower, curr, upper):
        A = numpy.matrix([[0] * (self.numx + 1)] * (self.numx + 1), dtype = numpy.float64)
        for i in range(self.numx + 1):
            a = self.alpha(i, time)
            b = self.beta(i, time)
            if i - 1 >= 0:
                A.itemset((i, i - 1), lower(a, b))
            A.itemset((i, i), curr(a, b))
            if i + 1 <= self.numx:
                A.itemset((i, i + 1), upper(a, b))
        return A

    def B_matrix(self, time, lower, curr, upper):
        A = numpy.matrix([[0] * (self.numx + 1)] * (self.numx + 1), dtype = numpy.float64)
        for i in range(self.numx + 1):
            a = self.alpha(i, time)
            b = self.beta(i, time)
            if i - 1 >= 0:
                A.itemset((i, i - 1), lower(a, b))
            A.itemset((i, i), curr(a, b))
            if i + 1 <= self.numx:
                A.itemset((i, i + 1), upper(a, b))
        return A

    # Stub methods
    def alpha(self, height, time):
        return 0

    def beta(self, height, time):
        return 0

    def a(self, t):
        return 0

    def b(self, t):
        return 0

    def initial_value_at_height(self, j):
        return 0

    def initial_value_at_top(self, i):
        return 0

    def initial_value_at_bottom(self, i):
        return 0

    def xi(self, s, t):
        return 0
