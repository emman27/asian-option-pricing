from floating_call import FloatingCall
import numpy
import math

class AsianDLFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, old_average, t0)
        self.xi_initial = self.xi(self.s0, 0)
        self.maxx = self.xi_initial * 3
        self.dx = self.maxx * 2 / self.numx
        self.j0 = round((self.xi_initial + self.maxx) / self.dx)
        self.set_boundary_conditions()

    def a(self, t):
        return t/(self.t0 + self.maxt) - self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(r*(self.maxt - t))

    def set_initial_boundary(self):
        for row in range(self.numx):
            self.grid.itemset((row, 0), self.initial_value_at_height(row))

    def initial_value_at_height(self, row):
        return max(self.t0 / (self.t0 + self.maxt) - self.maxx + row * self.dx, 0)

    def set_top_boundary(self):
        for col in range(self.numt):
            self.grid.itemset((self.numx - 1, col), self.initial_value_at_time(col))

    def initial_value_at_time(self, col):
        return 1 - 1/(self.r * (self.t0 + self.maxt)) + math.exp(-self.r * col * self.dt) * (
            1/(self.r * (self.t0 + self.maxt)) - (self.maxt - col * self.dt)/(self.t0 + self.maxt) + self.maxx
        )

    def alpha(self, height, time):
        return .5 * self.sigma**2 * (
            (self.maxt - (time - 0.5) * self.dt) / (self.t0 + self.maxt) +
            (self.numx/2 - height) * self.dx
        )**2 / (2 * self.dx**2) * self.dt

    def beta(self, height, time):
        return self.r * (
            (self.maxt - (time - 0.5) * self.dt) / (self.t0 + self.maxt) +
            (self.numx/2 - height) * self.dx
        ) / (4 * self.dx) * self.dt

    def A_matrix(self, time):
        A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i, time)
            b = self.beta(i, time)
            if i - 1 >= 0:
                A.itemset((i, i-1), a - b)
            A.itemset((i, i), 1 - 2 * a)
            if i + 1 < numx:
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
            if i + 1 < numx:
                B.itemset((i, i + 1), a + b)
        return B

    def solve(self):
        super().solve(lambda time: numpy.identity(self.numx) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.s0 * self.grid[self.j0, self.numt]

# Constants
numx = 200
numt = 400
maxt = 1
r = 0.1
s0 = 100
sigma = 0.3
#
# t0 = 0.1
# print('Expected: 9.85, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
# print('Expected: 9.34, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
# print('Expected: 8.84, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
#
# t0 = 0.3
# print('Expected: 10.70, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
# print('Expected: 9.05, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
# print('Expected: 7.61, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

t0 = 0.9
print('Expected: 10.38, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
print('Expected: 4.07, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
print('Expected: 1.03, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
