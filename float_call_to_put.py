import math

def put_price_floating(floating_call):
    return floating_call.solve() - floating_call.s0 + floating_call.s0 / (floating_call.r * (floating_call.t0 + floating_call.maxt)) * (1 - math.exp(-floating_call.r * floating_call.maxt))
