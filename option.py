from grid import Grid

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
