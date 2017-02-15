import math
import numpy

def montecarlo(curr, strike, r, sigma, maxt, numt, times):
  orig = curr
  dt = maxt / numt
  counts = []
  for j in range(times):
    curr = orig
    for i in range(numt):
      curr = curr * (1 + r * dt + sigma * math.sqrt(dt) * numpy.random.normal(0, 1))
    counts.append(curr - orig)
  return sum(counts) / len(counts)

print(montecarlo(100, 100, 0.09, 0.3, 1, 400, 100000))
