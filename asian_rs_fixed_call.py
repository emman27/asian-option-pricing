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
  def alpha(height):
    return .25 * self.sigma**2 * height**2 * self.dt

  def beta(height):
    return (height * self.r * self.dx + 1 / self.maxt) * (self.deft / (4 * self.dx))

  # The relation between two columns are given by
  # R = (I-B)^-1 * A * L
  # where they follow the original relation
  # R = AL + BR, with
  # A = [[1,0,0,0...]...[i-1th, ith, i+1th, 0, 0...], [0, i-1th, ith, i+1th, 0...]...[...0,0,0,1]]
  def A_matrix():
    A = numpy.matrix([[0] * numt] * numx, dtype = numpy.float64)
    for i in range(1, numx - 1):
      a = alpha(i)
      b = beta(i)
      A.itemset((i, i - 1), (a + b) * dt)
      A.itemset((i, i), (1/dt - 2*a) * dt)
      A.itemset((i, i + 1), (a - b) * dt)
    return A

  def B_matrix():
    B = numpy.matrix([[0] * numt] * numx, dtype = numpy.float64)
    B.itemset((0,0), 1)
    B.itemset((numx - 1, numt - 1), 1)
    for i in range(1, numx - 1):
      a = alpha(i)
      b = beta(i)
      A.itemset((i, i - 1), (a + b) * dt)
      A.itemset((i, i), (-2 * a) * dt)
      A.itemset((i, i + 1), (a - b) * dt)
    return B

print(AsianRSFixedCall(3, 1, 300, 100, 0.02, 0.3).grid)
