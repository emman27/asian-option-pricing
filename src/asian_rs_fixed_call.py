'''
Implements the Asian Option with the Rogers-Shi PDE
'''
import math
from src.option import Option

class AsianRSFixedCall(Option):
    '''
    Implements the Asian Option with the Rogers-Shi PDE
    '''
    def __init__(self, maxx, maxt, numx, numt, r, sigma):
        super().__init__(maxx, maxt, numx, numt, r, sigma)

        self.set_boundary_conditions()

    def set_boundary_conditions(self):
        '''
        Sets the boundary conditions for the grid based on PDE formulation
        '''
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

    def alpha(self, height):
        '''
        Calculates the alpha variable, abstraction of coefficients
        '''
        return .25 * self.sigma**2 * height**2 * self.dt

    def beta(self, height):
        '''
        Calculates the beta variable, abstraction of the coefficients
        '''
        return (height * self.r * self.dx + 1 / self.T) * (self.dt / (4 * self.dx))

    def mat_a(self, curr_col):
        '''
        Matrix A containing the coefficients applied to the current column
        '''
        pass

    def mat_b(self, curr_col):
        '''
        Matrix B containing the coefficients applied to the next column, aka the one we are solving
        '''
        pass

print(AsianRSFixedCall(3, 1, 300, 100, 0.02, 0.3).grid)
