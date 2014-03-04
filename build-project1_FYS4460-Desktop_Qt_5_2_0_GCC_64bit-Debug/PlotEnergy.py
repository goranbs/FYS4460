'''
                    PlotEnergy.py - project1_FYS4460

Should read output files from the MD program in FYS4460.
'''

import numpy as np

num_files = 499

filename = 'state%04g.txt' % 0      # filename = 'state0000.txt' initial file
file = open(filename,'r')

N = int(file.readline())            # number of atoms
secondline = file.readline()
somestring,t0 = secondline.split()
t0 = float(t0)                        # time in simulation
E_sys = []
time = []
time.append(t0) # first time value
E_sys.append(0) # first value

for i in range(num_files-1):
    i = i+1
    filename = 'state%04g.txt' % i      # filename = 'state0001.txt'
    file = open(filename,'r')
    
    N_atoms = int(file.readline())            # number of atoms
    secondline = file.readline()
    somestring,t = secondline.split()
    time.append(float(t))                        # time in simulation
    
    E_tot = 0
    j = 0
    for line_j in file:
        line = line_j.split()           # split line on whitespace
        E_tot += float(line[12])
        j += 1
        
    if j == N_atoms:
        print '%g OK' % i
    else:
        print 'Danger!! j != N'


    E_sys.append(E_tot)

##################################################
# plotting

import matplotlib.pyplot as plt

plt.figure()
plt.plot(time,E_sys)
plt.title('Energy of system')
plt.xlabel('unitless time')
plt.ylabel('unitless energy')
plt.legend('E(t)')
plt.show(True)
