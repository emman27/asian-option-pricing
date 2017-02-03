from option import Option
import math
import numpy

class AsianRSFixedCall(Option):
    def __init__(self, maxt, numx, numt, r, sigma, initial_price, strike):
        super().__init__(maxt, numx, numt, r, sigma, initial_price, strike)
        self.j0 = round(numx/3)
        self.xi_initial = self.xi(self.s0, 0)
        self.dx = self.strike / self.s0 / self.j0
        self.set_boundary_conditions()

    def set_bottom_boundary(self):
        '''
        Sets the boundary values for the FDS grid when the height is 0
        '''
        for col in range(self.numt + 1):
            time = col * self.dt
            self.grid.itemset((0, col), self.initial_zero_value_at_time(time))

    def initial_zero_value_at_time(self, time):
        '''
        Calculates the zero-height value of the FDS grid at given time @time
        '''
        top = (1 - math.exp(-self.r * time))
        bottom = (self.r * self.maxt)
        return top / bottom

    def a(self, t):
        return -self.q(t) * self.b(t)

    def b(self, t):
        return -math.exp(self.r*(self.maxt - t))

    # Methods to calculate the abstracted alpha and beta variables. See report.
    def alpha(self, height):
        return .25 * self.sigma**2 * height**2 * self.dt

    def beta(self, height):
        return (height * self.r * self.dx + 1 / self.maxt) * self.dt  / (4 * self.dx)

    # The relation between two columns are given by
    # (I-B) * R = A * L
    # where they follow the original relation
    # R = AL + BR. See report.
    def A_matrix(self):
        A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i)
            b = self.beta(i)
            try:
                A.itemset((i, i - 1), (a + b))
                A.itemset((i, i), (1 - 2*a))
                A.itemset((i, i + 1), (a - b))
            except:
                pass
        return A

    def B_matrix(self):
        B = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
        for i in range(self.numx):
            a = self.alpha(i)
            b = self.beta(i)
            try:
                B.itemset((i, i - 1), (a + b))
                B.itemset((i, i), (-2 * a))
                B.itemset((i, i + 1), (a - b))
            except:
                pass
        return B

    def solve(self):
        super().solve(numpy.identity(self.numx) - self.B_matrix(), self.A_matrix())
        return self.s0 * self.grid[self.j0, self.numt]

print(AsianRSFixedCall(1, 200, 400, 0.09, 0.05, 100, 90).solve())
