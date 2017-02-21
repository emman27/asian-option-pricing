from floating_call import FloatingCall
import math
import numpy

class AsianNewVecerFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, old_average, t0)
        self.maxx = self.xi_initial * 3
        self.dx = self.maxx / self.numx
        self.j0 = round(self.xi_initial / self.dx)
        self.set_boundary_conditions()

    def xi(self, s, t):
        return s / self.avr(t)

    def initial_value_at_height(self, row):
        return max(row * self.dx - 1, 0)

    def initial_value_at_top(self, col):
        return self.maxx - 1

    def alpha(self, height, time):
        return .5 * self.sigma**2 * (height * self.dx)**2 * (
            (1 - math.exp(-self.r * (height + 0.5) * self.dt)) /
            (self.r * (self.t0 + self.maxt)) *
            height * self.dx -
            1
        ) ** 2 * self.dt / (2 * self.dx**2)

    def A_matrix(self, time):
        return super().A_matrix(time, lambda a, b: a, lambda a, b: 1 - 2 * a, lambda a, b: a)

    def B_matrix(self, time):
        return super().B_matrix(time, lambda a, b: a, lambda a, b: - 2 * a, lambda a, b: a)

    def solve(self):
        super().solve(lambda time: numpy.identity(self.numx + 1) - self.B_matrix(time), lambda time: self.A_matrix(time))
        return self.s0 * self.grid[self.j0, self.numt]

if __name__ == "__main__":
    # Constants
    numx = 200
    numt = 400
    maxt = 1
    r = 0.1
    s0 = 100
    sigma = 0.3

    t0 = 0.1
    print('Expected: 9.87, Actual: ' + str(AsianNewVecerFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
    print('Expected: 9.35, Actual: ' + str(AsianNewVecerFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
    print('Expected: 8.84, Actual: ' + str(AsianNewVecerFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))

    t0 = 0.7
    print('Expected: 11.21, Actual: ' + str(AsianNewVecerFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
    print('Expected: 6.86, Actual: ' + str(AsianNewVecerFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
    print('Expected: 3.85, Actual: ' + str(AsianNewVecerFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
