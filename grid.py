import numpy

class Grid:
  def __init__(self, numx, numt):
    '''
    Initializes the Grid used for the Finite Difference Scheme
    '''
    self.matrix = numpy.matrix([[0] * numt] * numx)

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
