'''
Main class for Asian Option Pricing
Created by: Emmanuel Goh
'''

# Library for doing all the matrix stuff
import numpy as np
# Library for parsing command line arguments
import sys

def main():
  '''
  Main method in the function
  Arguments to pass, in the following order:
  @maxx: Maximum value of the stock at any point of time
  @maxt: Time to maturity of the stock, assuming current time t=0
  @numx: Number of nonzero price points of the stock
  @numt: Number of nonzero time points of monitoring
  Note that minimum price of a stock is simply 0
  '''
  maxx, maxt, numx, numt = sys.argv[1:] # First argument of argv is the script name


if __name__ == "__main__":
  main()
