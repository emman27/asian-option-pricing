from grid import Grid
from numpy.linalg import inv

class Option:
  def __init__(self, maxx, maxt, numx, numt, r, sigma):
    self.maxx = maxx
    self.T = maxt
    self.numx = numx
    self.numt = numt
    self.r = r
    self.sigma = sigma

    self.dx = maxx / float(numx)
    self.dt = maxt / float(numt)

    self.grid = Grid(maxx, maxt, numx, numt)

  def get_next_col(self, curr_col):
    '''
    Governed by the rule u_{i+1} = (I-B_i)^{-1} A_i u_i
    '''
    self.grid.set_col(inv(numpy.eye(self.grid.num_rows(), dtype=float) - B(curr_col)) * A(curr_col))

  def A(col):
    '''Abstract method'''
    pass

  def B(col):
    '''Abstract method'''
    pass
