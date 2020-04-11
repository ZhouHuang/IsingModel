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
    _MCTIME = 4       # total simulation times
    _MU = 1          # magnetic moment of the particle, borh magneton
    _H = 0           # magnetic field strength
    _T = 100         # system temperature, unit Kelvin
    _ENERGY = 0      # system energy
    _ENERGYVAR = 0   # system energy variance
    _MAGINTENSITY = 0      # system magnetic intensity
    _MAGINTENSITYVAR = 0   # system magnetic intensity variance
    _KBOLTZMANN = 1  # boltzmann constant
    
    __LABEL = {-1:'O', 1:'X'} # visulize spin
    
    def __init__(self,*args):
        if len(args) >= 2:
            print('initialing type ERROR! default type is all alined')
            print('USAGE:')
            return None
        elif len(args) == 0:
            self._STATE = mat([[2*np.pi for i in range(self._RAW)] for ii in range(self._COLUMN)])
            #_STATE = 2*random.randint(2, size=(N,N))-1 # random spin 
        elif len(args) == 1:
            print('args == 1')
        
    def calculateTotalEnergy(self):
        pass
    
    def getTotalEnergy(self):
        return self._ENERGY
    
    def getTotalEnergyVariance(self):
        return self._ENERGYVAR
    
    def calculateAdjacentEnergy(self):
        pass

    def calculateFieldEnergy(self):
        pass
    
    def calculateMagneticIntensity(self):
        pass
    
    def getMagneticIntensity(self):
        return self._MAGINTENSITY
    
    def getMagneticIntensityVariance(self):
        return self._MAGINTENSITYVAR
    
    def _randomFlip(self,*angle): # rotate angle, default is pi, unit rad
        raw = random.randint(self._RAW)
        column = random.randint(self._COLUMN)
        print(raw,column,angle)
        self._STATE[raw,column] = (self._STATE[raw,column] + np.pi) % (2*np.pi) # up-down flip
        #self._STATE[raw,column] = (self._STATE[raw,column] + angle) % (2*np.pi) # certain angle flip
        return (raw,column)

    def visulizeSpin(self):
        VS = mat([[self.__LABEL[np.cos(self._STATE[r,c])] \
                           for r in range(self._RAW) ] for c in range(self._COLUMN)] )
        print(VS)
        
    def simulate(self):
        for i in range(self._MCTIME):
            self._randomFlip()
    
def main():
    ising = twoDIsing()
    
    ising.visulizeSpin()
    print('-'*15,'fliping','-'*15)
    ising.simulate()
    ising.visulizeSpin()
    
    def plot():
        pass
        
if __name__ == '__main__':
    main()