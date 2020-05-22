# -*- coding: utf-8 -*-
import numpy as np

def main():
    dis = 'up'
    t = 0.1
    # with open('../results/data_{}_{}.dat'.format(dis,t),'w') as f:
    #     xi=0
    #     while(xi<len(x)):
    #         f.write(str(x[xi]) +'    ' + str(energy[xi]) + '    ' + str(magnetTensity[xi]) + '\n')
    #         xi +=1
    energy =[]
    magnetTensity =[]
    with open('../results/data_{}_{}.dat'.format(dis,t)) as f:
        line = f.readline()
        while line:
            # print(line.split('    '),type(line))
            step = float(line.split('    ')[0])
            if step >51:
                energy.append(float(line.split('    ')[1]))
                magnetTensity.append(float(line.split('    ')[2].rstrip('\n')))
            # print(step,energy,magnetTensity)
            # print(magnetTensity)
            line = f.readline()
            
    energy = np.array(energy)
    delta = []
    tau = []
    
    while (len(delta) < 6):
        energy_mean = sum(energy)/len(energy)
        sigma = np.sqrt((sum(energy**2)/len(energy) - energy_mean**2)/(len(energy)-1))
        # sigma_origin = sum(energy)/len(energy)
        print('energy',energy)
        print('energy mean',energy_mean)
        print('energy sigma',sigma)
    
    
        delta.append(sigma)
        print(energy.shape,len(energy))
        energy = np.array( [0.5*(energy[2*i]+energy[2*i+1]) for i in range(int(len(energy)/2))])
        print('energy',energy)
        print(energy.shape,len(energy))
        
        print('binning level',len(delta))
        print('delta',delta)
        tau.append(0.5* ((delta[-1]/delta[0])**2 - 1))
        print('tau',tau)
            
if __name__ == '__main__':
    main()