from option import Option
import math

class FixedCall(Option):
    def __init__(self, maxt, numx, numt, r, sigma, s0, strike):
        super().__init__(maxt, numx, numt, r, sigma, s0)

    def xi(self, s, t):
        return self.a(t) + self.b(t) * self.avr(t) / s
