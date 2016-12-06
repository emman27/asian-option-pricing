'''
Generic Option class, acts as a superclass for all Asian Options in this module
'''
import numpy
from numpy.linalg import inv
from grid import Grid

class Option:
    '''
    Generic Option class, acts as a superclass for all Asian Options in this module
    '''
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
        inverse = inv(numpy.eye(self.grid.num_rows(), dtype=float) - self.mat_b(curr_col))
        self.grid.set_col(inverse * self.mat_a(curr_col))

    def mat_a(self, col):
        '''Abstract method'''
        raise NotImplementedError('This method is an abstract method')

    def mat_b(self, col):
        '''Abstract method'''
        raise NotImplementedError('This method is an abstract method')
