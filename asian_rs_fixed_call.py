from grid import Grid
import math

class AsianRSFixedCall:
  def __init__(self, maxx, maxt, numx, numt, r, sigma):
    self.maxx = maxx
    self.T = maxt
    self.numx = numx
    self.numt = numt
    self.r = r
    self.sigma = sigma

    self.dx = maxx / numx
    self.dt = maxt / numt

    self.grid = Grid(maxx, maxt, numx, numt)

    self.set_boundary_conditions()

  def set_boundary_conditions(self):
    self.set_height_zero()

  def set_height_zero(self):
    '''
    Sets the boundary values for the FDS grid when the height is 0
    '''
    for col in range(self.numt):
      time = col * self.dt
      self.grid.set_value_at(0, time, self.initial_zero_value_at_time(time))

  def initial_zero_value_at_time(self, time):
    '''
    Calculates the zero-height value of the FDS grid at given time @time
    '''
    top = (1 - math.exp(-self.r * time))
    bottom = (self.r * self.T)
    return top / bottom
