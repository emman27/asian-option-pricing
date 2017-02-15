from .context import src
from src.grid import Grid
import unittest

import src.asian_rs_fixed_call
import src.asian_vecer_fixed_call
import src.asian_dl_fixed_call
import src.asian_new_vecer_fixed_call

import src.asian_rs_float_call
import src.asian_vecer_float_call
import src.asian_dl_float_call
import src.asian_new_vecer_float_call

class TestGridMethods(unittest.TestCase):
    def test_grid_size(self):
        '''
        Ensures that matrices are created with the correct dimensions
        '''
        grid = Grid(1, 1, 100, 99)
        self.assertEqual(100, grid.get_raw_matrix().shape[0])
        self.assertEqual(99, grid.get_raw_matrix().shape[1])

    def test_get_value(self):
        '''
        Ensures that the retrieval method of Grid works
        '''
        grid = Grid(1, 1, 50, 50)
        self.assertEqual(0, grid.get_value(0, 34))

    def test_set_value(self):
        '''
        Ensures that can correctly set a value. Dependent on get_value to work as well
        '''
        grid = Grid(1, 1, 10, 10)
        grid.set_value(3, 3, 22)
        self.assertEqual(22, grid.get_value(3, 3))

    def test_out_of_bounds_value(self):
        '''
        Ensures that an Exception is raised if the retrieval is invalid
        '''
        grid = Grid(1, 1, 1, 1)
        with self.assertRaises(IndexError):
            grid.get_value(1, 1)

    def test_get_column(self):
        grid = Grid(1, 1, 2, 2)
        grid.set_value(0, 0, 1)
        grid.set_value(1, 0, 1)
        self.assertEqual([[1], [1]], list(grid.get_col(0)))

    def test_set_column(self):
        pass

if __name__ == '__main__':
    unittest.main()
