# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:56:58 2020
@author: huangzhou
2D Ising model simulation
"""
import numpy as np
from numpy import random, mat, cos

class twoDIsing():
    # parameters definition
    _XMAX = 4     # X coordinate
    _YMAX = 4    # Y coordinate 
                  # note that X Y are on commen sense which equles the matrix element mat[x,y]
    _STATE = [[]]     # site spin, pi is up, pi/2 is down
    _J = -1          # J = -1 ferromagnetizem J = 1 antiferromagnetizem
    _MCTIME = 40       # total simulation times
    _MU = 1          # magnetic moment of the particle, borh magneton
    _H = 0           # magnetic field strength
    _T = 0.2         # system temperature, unit Kelvin
    _ENERGY = 0      # system energy
    _ENERGYVAR = 0   # system energy variance
    _FIELDENERGY = 0 # field energy
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
        self._ENERGYVAR = 0
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._ENERGY += self._calculateLocalEnergy(x,y)
        self._ENERGY = self._ENERGY/(self._XMAX * self._YMAX)
        #print('Total Energy = ',self._ENERGY)
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._ENERGYVAR += (self._calculateLocalEnergy(x,y) - self._ENERGY)**2
        self._ENERGYVAR = self._ENERGYVAR / (self._XMAX * self._YMAX)
    
    def getTotalEnergy(self):
        return self._ENERGY
    
    def getTotalEnergyVariance(self):
        return self._ENERGYVAR
    
    def _calculateLocalEnergy(self,x,y):
        # energy between adjacent sites
        # epsilon = sum{ sigma_xy * sigma_j }, j are neighboring
        x = x % self._XMAX
        y = y % self._YMAX
        top = [x, y - 1 if y>0 else self._YMAX - 1]
        bottom = [x, y + 1 if y < self._YMAX - 1 else 0]
        left = [x - 1 if x>0 else self._XMAX - 1 , y]
        right = [x + 1 if x < self._XMAX - 1 else 0 , y]
        epsilon = self._J * (cos(self._STATE[top[0],top[1]])*cos(self._STATE[x,y]) + \
                             cos(self._STATE[bottom[0],bottom[1]])*cos(self._STATE[x,y]) +\
                             cos(self._STATE[left[0],left[1]])*cos(self._STATE[x,y]) + \
                             cos(self._STATE[right[0],right[1]])*cos(self._STATE[x,y]) )
        epsilon = epsilon/2
        #print('neighboring energy = ',epsilon)
        return epsilon

    def _calculateFieldEnergy(self):
        self._FIELDENERGY = 0
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._FIELDENERGY += -1 * self._MU * self._H * cos(self._STATE[x,y])
        print('field energy = ',self._FIELDENERGY)
    
    def _calculateMagneticIntensity(self):
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._MAGINTENSITY += cos(self._STATE[x,y]) 
        self._MAGINTENSITY = self._MAGINTENSITY/(self._XMAX * self._YMAX)
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._MAGINTENSITYVAR += ( cos(self._STATE[x,y]) - self._MAGINTENSITY )**2
        self._MAGINTENSITYVAR = self._MAGINTENSITYVAR / (self._XMAX * self._YMAX)
    
    def getMagneticIntensity(self):
        return self._MAGINTENSITY
    
    def getMagneticIntensityVariance(self):
        return self._MAGINTENSITYVAR
    
    def _flip(self,x,y,*angle):
        #self._STATE[x,y] = (self._STATE[x,y] + np.pi) % (2*np.pi) # up-down flip
        self._STATE[x,y] = (self._STATE[x,y] + angle) % (2*np.pi) # certain angle flip


    def visulizeSpin(self):
        VS = mat([[self.__LABEL[cos(self._STATE[x,y])] \
                           for x in range(self._XMAX) ] for y in range(self._YMAX)] )
        print(VS)
        
    def simulate(self):
        for i in range(self._MCTIME):
            y = random.randint(self._YMAX)
            x = random.randint(self._XMAX)
            angle = np.pi
            
            self._calculateTotalEnergy()
            E_before = self._ENERGY
            self._flip(x,y,angle)
            self._calculateTotalEnergy()
            E_after = self._ENERGY
            print('E_before = ',E_before,'E_after = ',E_after)
            if E_after < E_before :
                print('-'*8,'site{} flipping'.format((x,y)))
            elif random.uniform(0.0,1.0) <= np.exp(-(E_after-E_before)/(self._KBOLTZMANN * self._T)):
                print('-'*8,'site{} flipping'.format((x,y)))
            else:
                print('-'*8,'site{} holds'.format((x,y)))
                self._flip(x,y,-angle)
    
def main():
    ising = twoDIsing()
    
    ising.visulizeSpin()
    ising.simulate()
    ising.visulizeSpin()
    
    def plot():
        pass
        
if __name__ == '__main__':
    main()