from option import Option
import math
import numpy

class FloatingCall(Option):
    def __init__(self, maxt, numx, numt, r, sigma, s0, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, s0, old_average, t0)

    def xi(self, s, t):
        return self.a(t) + self.b(t) * self.avr(t) / s


    def solve(self, left_multiplier, right_multiplier):
        for col in range(self.numt):
            l = self.grid[:, col]
            r = self.grid[:, col + 1]
            new = numpy.linalg.solve(left_multiplier(col), right_multiplier(col) * l)
            for row in range(self.numx):
                self.grid.itemset((row, col + 1), new[row])
