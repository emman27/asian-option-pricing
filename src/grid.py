import numpy

class Grid:
    def __init__(self, maxx, maxt, numx, numt):
        '''
        Initializes the Grid used for the Finite Difference Scheme
        '''
        self.matrix = numpy.matrix([[0] * numt] * numx, dtype=numpy.float64)
        self.maxx = maxx
        self.maxt = maxt

        # Set the difference between each point
        self.dx = maxx / float(numx)
        self.dt = maxt / float(numt)

    def get_raw_matrix(self):
        '''
        Returns the matrix raw by itself
        '''
        return self.matrix

    def num_rows(self):
        '''
        Returns the number of rows for the Grid
        '''
        return self.matrix.shape[0]

    def num_cols(self):
        '''
        Returns the number of columns for the Grid
        '''
        return self.matrix.shape[1]

    def set_value(self, row, col, val):
        '''
        Sets the value of the cell at (row, col) to val
        Input:
          @row: The row index to be changed. Note: Indexing starts from 0
          @col: The col index to be changed.
          @val: The new value for the cell
        '''
        self.matrix.itemset((row, col), val)

    def get_value(self, row, col):
        '''
        Returns the value of the cell at (row, cell)
        '''
        return self.matrix.item((row, col))

    def value_at(self, price, time):
      '''
        Returns the value found at the price and time given.
        Note: Floored
        '''
        return self.matrix.item((int(price / self.dx), int(time / self.dt)))

    def set_value_at(self, price, time, val):
        self.matrix.itemset((int(price / self.dx), int(time / self.dt)), val)

    def get_col(self, col_num):
        return self.matrix[:, col_num]

    def set_col(self, col_num, new_col):
       pass

    def __str__(self):
       return self.matrix.__str__()
