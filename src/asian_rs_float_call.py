from floating_call import FloatingCall
import numpy
import math

class AsianRSFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, old_average, t0)
        self.maxx = 2 # In this case, the generic abs(self.xi_initial * 3) is too small
        self.dx = self.maxx * 2 / self.numx
        self.j0 = round((self.xi_initial + self.maxx) / self.dx)
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
        super().solve(lambda time: b_mat, lambda time: a_mat)
        return self.s0 * self.grid[self.j0, self.numt]

# if __name__ == '__main__':
#     # Constants
#     numx = 200
#     numt = 400
#     maxt = 1
#     r = 0.1
#     s0 = 100
#     sigma = 0.3

#     t0 = 0.1
#     print('Expected: 9.87, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 9.34, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 8.85, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

#     t0 = 0.3
#     print('Expected: 10.72, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 9.06, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 7.61, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

#     t0 = 0.7
#     print('Expected: 11.19, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 6.86, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 3.85, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

#     t0 = 0.9
#     print('Expected: 10.40, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
#     print('Expected: 4.05, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
#     print('Expected: 1.03, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
