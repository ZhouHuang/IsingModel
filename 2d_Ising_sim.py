# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:56:58 2020
@author: huangzhou
2D Ising model simulation
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import random, mat, cos

class twoDIsing():
    # parameters definition
    _XMAX = 8     # X coordinate
    _YMAX = 8    # Y coordinate 
                  # note that X Y are on commen sense which equles the matrix element mat[x,y]
    _STATE = [[]]     # site spin, pi is up, pi/2 is down
    _J = -1          # J = -1 ferromagnetizem J = 1 antiferromagnetizem
    _MCTIME = 200       # total simulation times
    _MU = 1          # magnetic moment of the particle, borh magneton
    _H = 0           # magnetic field strength
    _T = 0.1         # system temperature, unit Kelvin
    _KBOLTZMANN = 1  # boltzmann constant
    _ENERGY = 0      # system energy
    _ENERGYVAR = 0   # system energy variance
    _FIELDENERGY = 0 # field energy
    _MAGINTENSITY = 0      # system magnetic intensity
    _MAGINTENSITYVAR = 0   # system magnetic intensity variance
    _ENERGYARRAY = []      # store the energy values
    _ENERGYVARARRAY = []
    _MAGINTENSITYARRAY = []
    _MAGINTENSITYVARARRAY = []
    
    __LABEL = {-1:'O', 1:'X'} # visulize spin
    
    def __init__(self,**kw):
        if len(kw) >= 2:
            print('initialing type ERROR! default type is all alined up')
            print('USAGE:[type=default,up,down,random]')
            return None
        elif len(kw) == 0:
            self._STATE = mat([[2*np.pi for i in range(self._YMAX)] for ii in range(self._XMAX)])
            self._T = 1
        else:
            for k in kw:
                if k=='type':
                    if kw[k]=='down':
                        self._STATE = mat([[np.pi for i in range(self._YMAX)] for ii in range(self._XMAX)])
                    elif kw[k]=='up' or kw[k]=='default':
                        self._STATE = mat([[2*np.pi for i in range(self._YMAX)] for ii in range(self._XMAX)])
                    elif kw[k]=='random':
                        self._STATE = mat( np.pi*random.randint(2, size=(self._XMAX,self._YMAX)))
                    else:
                        pass
    
    def getMCTime(self):
        return self._MCTIME
    def getSize(self):
        return (self._XMAX,self._YMAX)
        
    def _calculateTotalEnergy(self):
        self._ENERGY = 0
        self._ENERGYVAR = 0
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._ENERGY += self._calculateLocalEnergy(x,y)
        self._calculateFieldEnergy()
        self._ENERGY += self._FIELDENERGY
        self._ENERGY = self._ENERGY/(self._XMAX * self._YMAX)
        #print('Total Energy = ',self._ENERGY)
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._ENERGYVAR += (self._calculateLocalEnergy(x,y) - self._ENERGY)**2
        self._ENERGYVAR = self._ENERGYVAR / (self._XMAX * self._YMAX)
    
    def getTotalEnergy(self):
        return self._ENERGYARRAY
    
    def getTotalEnergyVariance(self):
        return self._ENERGYVARARRAY
    
    def _calculateLocalEnergy(self,x,y):
        # energy between adjacent sites
        # epsilon = sum{ sigma_xy * sigma_j }, j are neighboring
        epsilon = 0
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
        #print('field energy = ',self._FIELDENERGY)
    
    def _calculateMagneticIntensity(self):
        self._MAGINTENSITY = 0
        self._MAGINTENSITYVAR = 0
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._MAGINTENSITY += cos(self._STATE[x,y]) 
        self._MAGINTENSITY = self._MAGINTENSITY/(self._XMAX * self._YMAX)
        for x in range(self._XMAX):
            for y in range(self._YMAX):
                self._MAGINTENSITYVAR += ( cos(self._STATE[x,y]) - self._MAGINTENSITY )**2
        self._MAGINTENSITYVAR = self._MAGINTENSITYVAR / (self._XMAX * self._YMAX)
    
    def getMagneticIntensity(self):
        return self._MAGINTENSITYARRAY
    
    def getMagneticIntensityVariance(self):
        return self._MAGINTENSITYVARARRAY
    
    def _flip(self,x,y,*angle):
        #self._STATE[x,y] = (self._STATE[x,y] + np.pi) % (2*np.pi) # up-down flip
        self._STATE[x,y] = (self._STATE[x,y] + angle) % (2*np.pi) # certain angle flip


    def visulizeSpin(self):
        VS = mat([[self.__LABEL[cos(self._STATE[x,y])] \
                           for x in range(self._XMAX) ] for y in range(self._YMAX)] )
        print(VS)
        
    def simulate(self,**kw):
        self._calculateTotalEnergy()
        self._calculateMagneticIntensity()
        self._ENERGYARRAY.append(self._ENERGY)
        self._ENERGYVARARRAY.append(self._ENERGYVAR)
        self._MAGINTENSITYARRAY.append(self._MAGINTENSITY)
        self._MAGINTENSITYVARARRAY.append(self._MAGINTENSITYVAR)
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
                print('^'*8,'site{} flipping'.format((x,y)))
            elif random.uniform(0.0,1.0) <= np.exp(-(E_after-E_before)/(self._KBOLTZMANN * self._T)):
                print('^'*8,'site{} flipping'.format((x,y)))
            else:
                print('-'*8,'site{} holds'.format((x,y)))
                self._flip(x,y,-angle)
            self._calculateMagneticIntensity()
            self._ENERGYARRAY.append(self._ENERGY)
            self._ENERGYVARARRAY.append(self._ENERGYVAR)
            self._MAGINTENSITYARRAY.append(self._MAGINTENSITY)
            self._MAGINTENSITYVARARRAY.append(self._MAGINTENSITYVAR)
    
def main():
    ising = twoDIsing(type='down')
    
    ising.visulizeSpin()
    ising.simulate()
    ising.visulizeSpin()
    magnetTensity = ising.getMagneticIntensity()
    energy = ising.getTotalEnergy()
    magnetTensityVar = ising.getMagneticIntensityVariance()
    energyVar = ising.getTotalEnergyVariance()
    
    # smooth sample plot
    energy = energy[::int(0.01*len(energy))]
    energyVar = energyVar[::int(0.01*len(energyVar))]
    magnetTensity = magnetTensity[::int(0.01*len(magnetTensity))]
    magnetTensityVar = magnetTensityVar[::int(0.01*len(magnetTensityVar))]
    x = np.linspace(0,len(energy),len(energy))
    
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=False,figsize=[6,10])
    ax1.errorbar(x=x,y=energy,yerr=energyVar,fmt='o',ecolor='r',color='b')
    ax1.set_title('energy evolution')
    ax2.errorbar(x=x,y=magnetTensity,yerr=magnetTensityVar,fmt='o',ecolor='r',color='b')
    ax2.set_title('magnetic momentum evolution')

    
    '''#detail plot
    x = np.linspace(0, ising.getMCTime()+1,ising.getMCTime()+1)
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=False,figsize=[6,10])
    ax1.errorbar(x=x,y=energy,yerr=energyVar,fmt='o',ecolor='r',color='b')
    ax1.set_title('energy evolution')
    ax2.errorbar(x=x,y=magnetTensity,yerr=magnetTensityVar,fmt='o',ecolor='r',color='b')
    ax2.set_title('magnetic momentum evolution')
    '''


if __name__ == '__main__':
    main()