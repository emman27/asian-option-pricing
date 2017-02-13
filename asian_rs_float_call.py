from floating_call import FloatingCall
import numpy
import math

class AsianRSFloatCall(FloatingCall):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, old_average, t0)
        self.xi_initial = self.xi(self.s0, 0)
        self.maxx = abs(self.xi_initial * 3)
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

    def alpha(self, height):
        return 1/2 * self.sigma**2 * (-self.numx / 2 + height)**2 * self.dx**2 / (2*self.dx**2) * self.dt

    def beta(self, height):
        return -( self.r * (-self.numx / 2 + height) * self.dx + 1 / (self.t0 + self.maxt) ) / (4 * self.dx) * self.dt

    def A_matrix(self):
        A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i)
            b = self.beta(i)
            if i - 1 >= 0:
                A.itemset((i, i - 1), (a - b))
            A.itemset((i, i), (1 - 2*a))
            if i + 1 < self.numx:
                A.itemset((i, i + 1), (a + b))
        return A

    def B_matrix(self):
        B = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i)
            b = self.beta(i)
            if i - 1 >= 0:
                B.itemset((i, i - 1), (a - b))
            B.itemset((i, i), (-2 * a))
            if i + 1 < self.numx:
                B.itemset((i, i + 1), (a + b))
        return B

    def solve(self):
        a_mat = self.A_matrix()
        b_mat = numpy.identity(self.numx) - self.B_matrix()
        super().solve(lambda time: b_mat, lambda time: a_mat)
        return self.s0 * self.grid[self.j0, self.numt]

if __name__ == '__main__':
    # Constants
    numx = 200
    numt = 400
    maxt = 1
    r = 0.1
    s0 = 100
    sigma = 0.3

    t0 = 0.1
    print('Expected: 9.87, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
    # print('Expected: 9.34, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
    # print('Expected: 8.85, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
    #
    # t0 = 0.3
    # print('Expected: 10.72, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
    # print('Expected: 9.06, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
    # print('Expected: 7.61, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
    #
    # t0 = 0.7
    # print('Expected: 11.19, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
    # print('Expected: 6.86, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
    # print('Expected: 3.85, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
    #
    # t0 = 0.9
    # print('Expected: 10.40, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 90, t0).solve()))
    # print('Expected: 4.05, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 100, t0).solve()))
    # print('Expected: 1.03, Actual: ' + str(AsianRSFloatCall(maxt - t0, numx, numt, r, sigma, s0, 110, t0).solve()))
