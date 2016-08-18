'''
Main class for Asian Option Pricing
Created by: Emmanuel Goh
'''

# Library for doing all the matrix stuff
import numpy
# Library for parsing command line arguments
import sys

def main():
  '''
  Main method in the function
  Arguments to pass, in the following order:
  @numx: Number of different price points of the stock
  @numt: Number of different time points of monitoring
  @maxx: Maximum price of the stock at any point of time
  @maxt: Time to maturity of the stock, assuming current time t=0
  '''
  print sys.argv

if __name__ == "__main__":
  main()
