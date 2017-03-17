from floating_call import FloatingCall
import numpy
import math

class AsianDLFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx * 2, numt, r, sigma, initial_price, old_average, t0)
        self.dx = self.maxx * 2 / self.numx
        self.j0 = int((self.xi_initial + self.maxx) / self.dx)
        self.xi_min = -self.maxx
        self.set_boundary_conditions()

    def a(self, t):
        return t/(self.t0 + self.maxt) - self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    def initial_value_at_height(self, row):
        return max(self.t0 / (self.t0 + self.maxt) - self.maxx + row * self.dx, 0)

    def initial_value_at_top(self, col):
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
        return super().A_matrix(time, lambda a, b: a - b, lambda a, b: 1 - 2 * a, lambda a, b: a + b)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a - b, lambda a, b: - 2 * a, lambda a, b: a + b)

    def solve(self):
        return super().solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))

# if __name__ == "__main__":
#     # Constants
#     numx = 200
#     numt = 400
#     maxt = 1
#     r = 0.1
#     s0 = 100
#     sigma = 0.3

#     t0 = 0.1
#     print('Expected: 9.85, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 9.34, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 8.84, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

#     t0 = 0.3
#     print('Expected: 10.70, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 9.05, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 7.61, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

#     t0 = 0.9
#     print('Expected: 10.38, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 4.07, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 1.03, Actual: ' + str(AsianDLFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
