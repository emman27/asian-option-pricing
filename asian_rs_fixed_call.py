from option import Option
import math
import numpy

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
  def alpha(self, height):
    return .25 * self.sigma**2 * height**2

  def beta(self, height):
    return (height * self.r * self.dx + 1 / self.T)  / (4 * self.dx)

  # The relation between two columns are given by
  # R = (I-B)^-1 * A * L
  # where they follow the original relation
  # R = AL + BR, with
  # A = [[1,0,0,0...]...[i-1th, ith, i+1th, 0, 0...], [0, i-1th, ith, i+1th, 0...]...[...0,0,0,1]]
  def A_matrix(self):
    A = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
    for i in range(1, self.numx - 1):
      a = self.alpha(i)
      b = self.beta(i)
      A.itemset((i, i - 1), (a + b) * self.dt)
      A.itemset((i, i), (1/self.dt - 2*a) * self.dt)
      A.itemset((i, i + 1), (a - b) * self.dt)
    return A

  def B_matrix(self):
    B = numpy.matrix([[0] * self.numx] * self.numx, dtype = numpy.float64)
    B.itemset((0,0), 1)
    B.itemset((self.numx - 1, self.numx - 1), 1)
    for i in range(1, self.numx - 1):
      a = self.alpha(i)
      b = self.beta(i)
      B.itemset((i, i - 1), (a + b) * self.dt)
      B.itemset((i, i), (-2 * a) * self.dt)
      B.itemset((i, i + 1), (a - b) * self.dt)
    return B

  def solve(self):
    A = self.A_matrix()
    B = self.B_matrix()

    top_coeff = B[1, 0]
    btm_coeff = B[-2, -1]

    a_mat = self.A_matrix()[1:self.numx-1, 1:self.numx-1]
    b_mat = self.B_matrix()[1:self.numx-1, 1:self.numx-1]
    for i in range(self.numt - 1):
      L = self.grid.get_raw_matrix()[1:self.numx-1, i]
      k = numpy.matrix([[0]] * (self.numx - 2), dtype=numpy.float64)
      k.itemset((0,0), self.grid.get_raw_matrix()[0, i+1] * top_coeff)
      k.itemset((-1, 0), self.grid.get_raw_matrix()[self.numx - 1, i+1] * btm_coeff)
      new = (numpy.identity(self.numx-2) - b_mat).getI() * (a_mat * L + k)
      for j in range(self.numx - 2):
        self.grid.set_value(j + 1, i + 1, new[j, 0])
    return self.grid.get_raw_matrix()

print(AsianRSFixedCall(3, 1, 30, 10, 0.02, 0.3).solve())
