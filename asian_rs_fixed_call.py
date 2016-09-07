from option import Option
import math

class AsianRSFixedCall(Option):
  def __init__(self, maxx, maxt, numx, numt, r, sigma):
    super().__init__(maxx, maxt, numx, numt, r, sigma)

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


  # Methods to calculate the abstracted alpha and beta variables.
  # Math needs to be confirmed
  def self.alpha(height):
    return .25 * self.sigma**2 * height**2 * self.dt

  def self.beta(height):
    return (height * self.r * self.dx + 1 / self.maxt) * (self.dt / (4 * self.dx))

print(AsianRSFixedCall(3, 1, 300, 100, 0.02, 0.3).grid)
