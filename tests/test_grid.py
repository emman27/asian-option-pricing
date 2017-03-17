import unittest
import math

from .context import src

from src.option import Option

from src.fixed_call import FixedCall
from src.floating_call import FloatingCall

from src.asian_rs_fixed_call import AsianRSFixedCall
from src.asian_vecer_fixed_call import AsianVecerFixedCall
from src.asian_dl_fixed_call import AsianDLFixedCall
from src.asian_new_fixed_call import AsianNewFixedCall
from src.asian_new_vecer_fixed_call import AsianNewVecerFixedCall

from src.asian_rs_float_call import AsianRSFloatCall
from src.asian_vecer_float_call import AsianVecerFloatCall
from src.asian_dl_float_call import AsianDLFloatCall
from src.asian_new_float_call import AsianNewFloatCall
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

fix = FixedCall(MAXT, NUMX, NUMT, R, SIGMAS[0], S0, STRIKES[0])
flt = FloatingCall(MAXT, NUMX, NUMT, R, SIGMAS[0], S0, CURRENT_AVERAGES[0], T0S[0])

class TestOptions(unittest.TestCase):

    def test_max_xi_and_dx_positive(self):
        # Options to test
        fixed = [AsianRSFixedCall, AsianVecerFixedCall, AsianDLFixedCall, AsianNewFixedCall, AsianNewVecerFixedCall]
        floating = [AsianRSFloatCall, AsianVecerFloatCall, AsianDLFloatCall, AsianNewFloatCall, AsianNewVecerFloatCall]

        for option in fixed:
            for s in SIGMAS:
                for k in STRIKES:
                    opt = option(MAXT, NUMX, NUMT, R, s, S0, k)
                    self.assertGreaterEqual(opt.maxx, 0, msg=str(option))
                    self.assertGreaterEqual(opt.dx, 0)
        for option in floating:
            for s in SIGMAS_FLOAT:
                for t0 in T0S:
                    for avr in CURRENT_AVERAGES:
                        opt = option(MAXT - t0, NUMX, NUMT, R_FLOAT, s, S0, avr, t0)
                        self.assertGreaterEqual(opt.maxx, 0, msg=str(option))
                        self.assertGreaterEqual(opt.dx, 0)

    def test_xi_initial_and_j0(self):
        # Options to test. No restrictions on other options
        fixed = [AsianRSFixedCall, AsianDLFixedCall, AsianNewVecerFixedCall]
        floating = [AsianNewVecerFloatCall]
        the_rest_fixed = [AsianVecerFixedCall, AsianNewFixedCall]
        the_rest_floating = [AsianRSFloatCall, AsianDLFloatCall, AsianNewVecerFloatCall]

        for option in fixed:
            for s in SIGMAS:
                for k in STRIKES:
                    call = option(MAXT, NUMX, NUMT, R, s, S0, k)
                    self.assertGreaterEqual(call.xi_initial, 0, msg=str(option))
                    self.assertLess(abs(round(call.xi_initial, 4) - (round(call.j0 * call.dx,4))), call.dx, msg=str(call) + ': ' + str(call.xi_initial) + ', ' + str(call.dx * call.j0) + ', ' + str(call.j0))
        self.assertGreaterEqual(0, AsianNewFixedCall(MAXT, NUMX, NUMT, R, s, S0, k).xi_initial, msg=str(option))
        for option in floating:
            for s in SIGMAS_FLOAT:
                for t0 in T0S:
                    for avr in CURRENT_AVERAGES:
                        call = option(MAXT - t0, NUMX, NUMT, R_FLOAT, s, S0, avr, t0)
                        self.assertGreaterEqual(call.xi_initial, 0, msg=str(option))
                        self.assertLess(call.xi_initial - call.dx * call.j0, call.dx)
        for option in the_rest_fixed:
            for s in SIGMAS:
                for k in STRIKES:
                    call = option(MAXT, NUMX, NUMT, R, s, S0, k)
                    self.assertGreaterEqual(call.xi_initial, -call.maxx, msg=str(option))
                    self.assertGreaterEqual(call.maxx, call.xi_initial, msg=str(option))
                    self.assertLess(abs(round(call.xi_initial, 4) - (round(call.j0 * call.dx,4) - call.maxx)), call.dx, msg=str(call) + ': ' + str(call.xi_initial) + ', ' + str(call.dx * call.j0) + ', ' + str(call.j0))
        for option in the_rest_floating:
            for s in SIGMAS_FLOAT:
                for t0 in T0S:
                    for avr in CURRENT_AVERAGES:
                        self.assertGreaterEqual(call.xi_initial, -call.maxx, msg=str(option))
                        self.assertGreaterEqual(call.maxx, call.xi_initial, msg=str(option))
                        self.assertLess(abs(round(call.xi_initial, 4) - (round(call.j0 * call.dx,4) - call.maxx)), call.dx, msg=str(call) + ': ' + str(call.xi_initial) + ', ' + str(call.dx * call.j0) + ', ' + str(call.j0))


    def test_initializers(self):
        fixed = FixedCall(MAXT, NUMX, NUMT, R, .3, S0, 100)
        self.assertEqual(fixed.t0, 0)
        self.assertEqual(fixed.old_average, 0)

        fl = FloatingCall(MAXT - 0.1, NUMX, NUMT, R, .3, S0, 100, 0.1)
        self.assertEqual(fl.t0, 0.1)
        self.assertEqual(fl.old_average, 100)
        self.assertEqual(fl.maxt, 0.9)

    def test_tridiagonal(self):
        CONST_A = 1
        CONST_B = 2
        class TestOption(Option):
            def __init__(self, maxt, numx, numt, r, sigma, s0, avr, t0):
                super().__init__(maxt, numx, numt, r, sigma, s0, avr, t0)
            def alpha(self, i, j):
                return CONST_A
            def beta(self, i, j):
                return CONST_B
        opt = TestOption(MAXT, NUMX, NUMT, R, 0.3, S0, 100, 0.1)
        a = opt.A_matrix(1, lambda a, b: a, lambda a, b: b, lambda a, b: a + b)
        b = opt.B_matrix(1, lambda a, b: a, lambda a, b: b, lambda a, b: a + b)
        for mat in (a, b):
            for i in range(len(mat)):
                for j in range(mat[i, :].size):
                    if i - 1 == j:
                        self.assertEqual(mat[i, j], CONST_A)
                    elif i == j:
                        self.assertEqual(mat[i, j], CONST_B)
                    elif i + 1 == j:
                        self.assertEqual(mat[i, j], CONST_A + CONST_B)
                    else:
                        self.assertEqual(mat[i, j], 0)

    def test_q(self):
        self.assertEqual(fix.q(0), (1 - math.exp(-R * (MAXT - 0))) / (R * MAXT))
        self.assertEqual(flt.q(0), (1 - math.exp(-R * (MAXT - 0))) / (R * (T0S[0] + MAXT)))

    def test_average(self):
        self.assertEqual(fix.avr(0), fix.q(0) * S0 + math.exp(-R * MAXT) * 0)
        self.assertEqual(flt.avr(0), flt.q(0) * S0 + math.exp(-R * MAXT) * T0S[0] * CURRENT_AVERAGES[0] / (T0S[0] + MAXT))

    def test_stubs(self):
        opt = Option(MAXT, NUMX, NUMT, R, 0.3, S0, 100, 0.1)
        self.assertEqual(opt.alpha(3, 7), 0)
        self.assertEqual(opt.beta(3, 7), 0)
        self.assertEqual(opt.a(2), 0)
        self.assertEqual(opt.b(3), 0)
        self.assertEqual(opt.initial_value_at_top(124), 0)
        self.assertEqual(opt.initial_value_at_bottom(124), 0)
        self.assertEqual(opt.initial_value_at_height(124), 0)
        self.assertEqual(opt.xi(124, 32), 0)

if __name__ == '__main__':
    unittest.main()
