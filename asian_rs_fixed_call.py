from grid import Grid
import math

class AsianRSFixedCall:
  def __init__(self, maxx, maxt, numx, numt, r):
    self.maxx = maxx
    self.T = maxt
    self.numx = numx
    self.numt = numt
    self.r = r

    self.dx = maxx / numx
    self.dt = maxt / numt

    self.grid = Grid(maxx, maxt, numx, numt)

    self.set_boundary_conditions()

  def set_boundary_conditions(self):
    self.set_height_zero()

  def set_height_zero(self):
    for col in range(self.numt):
      time = col * self.dt
      self.grid.set_value_at(0, time, self.initial_zero_value_at_time(time))

  def initial_zero_value_at_time(self, time):
    top = (1 - math.exp(-self.r * time))
    bottom = (self.r * self.T)
    return top / bottom
