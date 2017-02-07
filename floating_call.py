from option import Option
import math
import numpy

class FloatingCall(Option):
    def __init__(self, maxt, numx, numt, r, sigma, s0, old_average, t0):
        super().__init__(maxt, numx, numt, r, sigma, s0, old_average, t0)

    def xi(self, s, t):
        return self.a(t) + self.b(t) * self.avr(t) / s
