from option import Option
import math

class FixedCall(Option):
    def __init__(self, maxt, numx, numt, r, sigma, s0, strike):
        super().__init__(maxt, numx, numt, r, sigma, s0, 0, 0)
        # FOr fixed calls, begin previous average and origin time at 0 (why?)
        self.strike = strike

    def xi(self, s, t):
        return self.a(t) + self.b(t) * (self.avr(t) - self.strike * math.exp(-self.r * (self.maxt - t))) / s
