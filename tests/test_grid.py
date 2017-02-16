from .context import src
import unittest

from src.asian_rs_fixed_call import AsianRSFixedCall
from src.asian_vecer_fixed_call import AsianVecerFixedCall
from src.asian_dl_fixed_call import AsianDLFixedCall
from src.asian_new_vecer_fixed_call import AsianNewVecerFixedCall

from src.asian_rs_float_call import AsianRSFloatCall
from src.asian_vecer_float_call import AsianVecerFloatCall
from src.asian_dl_float_call import AsianDLFloatCall
from src.asian_new_vecer_float_call import AsianNewVecerFloatCall


NUMX = 400
NUMT = 400
MAXT = 1
R = 0.09
S0 = 100
SIGMAS = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
STRIKES = [90, 95, 100, 105, 110]

R_FLOAT = 0.1
S0_FLOAT = 100
SIGMAS_FLOAT = [0.3, 0.5]
T0S = [0.1, 0.3, 0.5, 0.7, 0.9]
CURRENT_AVERAGES = [90, 100, 110]

class TestOptions(unittest.TestCase):

    def test_max_xi_positive(self):
        # Options to test
        fixed = [AsianRSFixedCall, AsianVecerFixedCall, AsianDLFixedCall, AsianNewVecerFixedCall]
        floating = [AsianRSFloatCall, AsianVecerFloatCall, AsianDLFloatCall, AsianNewVecerFloatCall]

        for option in fixed:
            for s in SIGMAS:
                for k in STRIKES:
                    self.assertGreaterEqual(option(MAXT, NUMX, NUMT, R, s, S0, k).maxx, 0, msg=str(option))
        for option in floating:
            for s in SIGMAS_FLOAT:
                for t0 in T0S:
                    for avr in CURRENT_AVERAGES:
                        self.assertGreaterEqual(option(MAXT - t0, NUMX, NUMT, R_FLOAT, s, S0, avr, t0).maxx, 0, msg=str(option))

    def test_xi_initial_within_range(self):
        # Options to test
        fixed = [AsianRSFixedCall, AsianDLFixedCall, AsianNewVecerFixedCall]
        floating = [AsianNewVecerFloatCall]

        for option in fixed:
            for s in SIGMAS:
                for k in STRIKES:
                    self.assertGreaterEqual(option(MAXT, NUMX, NUMT, R, s, S0, k).xi_initial, 0, msg=str(option))
        for option in floating:
            for s in SIGMAS_FLOAT:
                for t0 in T0S:
                    for avr in CURRENT_AVERAGES:
                        self.assertGreaterEqual(option(MAXT - t0, NUMX, NUMT, R_FLOAT, s, S0, avr, t0).xi_initial, 0, msg=str(option))


if __name__ == '__main__':
    unittest.main()
