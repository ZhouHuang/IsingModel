# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:56:58 2020
@author: huangzhou
2D Ising model simulation
"""
import numpy as np
from numpy import random, mat, array

class twoDIsing():
    # parameters definition
    _RAW = 10
    _COLUMN = 10     # the size of the interaction plane
    _STATE = [[]]     # site spin, pi is up, pi/2 is down
    _J = -1          # J = -1 ferromagnetizem J = 1 antiferromagnetizem
    _MCS = 100       # total simulation times
    _MU = 1          # magnetic moment of the particle, borh magneton
    _H = 0           # magnetic field strength
    _T = 100         # system temperature, unit Kelvin
    _ENERGY = 0      # system energy
    _ENERGYVAR = 0   # system energy variance
    _MAGINTENSITY = 0      # system magnetic intensity
    _MAGINTENSITYVAR = 0   # system magnetic intensity variance
    _KBOLTZMANN = 1  # boltzmann constant
    
    def __init__(self,*args):
        if len(args) >= 2:
            print('initialing type ERROR! default type is all alined')
            print('USAGE:')
            return None
        elif len(args) == 0:
            self._STATE = mat([[np.pi for i in range(self._RAW)] for ii in range(self._COLUMN)])
            #_STATE = 2*np.random.randint(2, size=(N,N))-1 # random spin 
        elif len(args) == 1:
            print('args == 1')
        
    def calculateTotalEnergy():
        pass
    def calculateAdjacentEnergy():
        pass
    def calculateFieldEnergy():
        pass
    def calculateMagneticIntensity():
        pass
    def randomFlip():
        pass
    def simulate():
        pass
    
def main():
    ising = twoDIsing()
    print(ising._STATE)
    
    def plot():
        pass
        
if __name__ == '__main__':
    main()