# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:56:58 2020
@author: huangzhou
2D Ising model simulation
"""
import numpy as np
from numpy import random, mat, array, cos, sin

class twoDIsing():
    # parameters definition
    _XMAX = 3     # X coordinate
    _YMAX = 10    # Y coordinate 
                  # note that X Y are on commen sense which equles the matrix element mat[x,y]
    _STATE = [[]]     # site spin, pi is up, pi/2 is down
    _J = -1          # J = -1 ferromagnetizem J = 1 antiferromagnetizem
    _MCTIME = 2       # total simulation times
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
            self._STATE = mat([[2*np.pi for i in range(self._YMAX)] for ii in range(self._XMAX)])
            #_STATE = 2*random.randint(2, size=(N,N))-1 # random spin 
        elif len(args) == 1:
            print('args == 1')
        
    def _calculateTotalEnergy(self):
        self._ENERGY = 0
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                energy += self._calculateAdjacentEnergy(x,y)
        self._ENERGY = self._ENERGY/2
        print('Total Energy = ',self._ENERGY)
        return self._ENERGY
    
    def getTotalEnergy(self):
        return self._ENERGY
    
    def getTotalEnergyVariance(self):
        return self._ENERGYVAR
    
    def _calculateAdjacentEnergy(self,x,y):
        x = x % self._XMAX
        y = y % self._YMAX
        top = [x, y - 1 if y>0 else self._YMAX - 1]
        bottom = [x, y + 1 if y < self._YMAX - 1 else 0]
        left = [x - 1 if x>0 else self._XMAX - 1 , y]
        right = [x + 1 if x < self._XMAX - 1 else 0 , y]
        # energy between adjacent sites
        #epsilon = sum{ sigma_xy * sigma_j }, j are neighboring
        epsilon = self._J * (cos(self._STATE[top[0],top[1]])*cos(self._STATE[x,y]) + \
                             cos(self._STATE[bottom[0],bottom[1]])*cos(self._STATE[x,y]) +\
                             cos(self._STATE[left[0],left[1]])*cos(self._STATE[x,y]) + \
                             cos(self._STATE[right[0],right[1]])*cos(self._STATE[x,y]) )
        print('neighboring energy = ',epsilon)
        return epsilon


    def _calculateFieldEnergy(self):
        pass
    
    def _calculateMagneticIntensity(self):
        pass
    
    def getMagneticIntensity(self):
        return self._MAGINTENSITY
    
    def getMagneticIntensityVariance(self):
        return self._MAGINTENSITYVAR
    
    def _randomFlip(self,*angle): # rotate angle, default is pi, unit rad
        y = random.randint(self._YMAX)
        x = random.randint(self._XMAX)
        print('site{} flipping{}'.format((x,y),angle))
        #self._STATE[x,y] = (self._STATE[x,y] + np.pi) % (2*np.pi) # up-down flip
        self._STATE[x,y] = (self._STATE[x,y] + angle) % (2*np.pi) # certain angle flip
        return (x,y)

    def visulizeSpin(self):
        VS = mat([[self.__LABEL[cos(self._STATE[x,y])] \
                           for x in range(self._XMAX) ] for y in range(self._YMAX)] )
        print(VS)
        
    def simulate(self):
        for i in range(self._MCTIME):
            self._randomFlip(np.pi)
    
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