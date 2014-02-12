'''
This program should read velocity data from stateXXX.txt file generated
by the MD program main.cpp and plot histograms of the velocity distribution
of the particles.
The particles should eventually after some time have a velocity distribution
similar to a Maxwell-Bolzmann distribution.
'''
import numpy as np
import os, sys


filename = 'state%03g.txt' % 0      # filename = 'state000.txt' initial file
file = open(filename,'r')

N = int(file.readline())            # number of atoms
secondline = file.readline()
somestring,t = secondline.split()
t = float(t)                        # time in simulation

files = []
for i in range(199-1):
    i = i+1
    filename = 'state%03g.txt' % i      # filename = 'state000.txt'
    file = open(filename,'r')
    
    N_atoms = int(file.readline())            # number of atoms
    secondline = file.readline()
    somestring,t = secondline.split()
    t = float(t)                        # time in simulation
    
    v_x = [] # speed in x -direction
    v_y = [] #          y -direction
    v_z = [] #          z -direction
    j = 0
    for line_j in file:
        line = line_j.split()           # split line on whitespace
        v_x.append(float(line[4]))
        v_y.append(float(line[5]))
        v_z.append(float(line[6]))
        j += 1
        
    if j == N_atoms:
        print '%g OK' % i
    else:
        print ' j != N, j=%g, N=%g' % j,N
    
    # plotting    
    import matplotlib.pyplot as plt
    bins = 10
    
    fig = plt.figure(figsize=(5,5))
    plt.subplot(3,1,1)
    plt.title('velocity distribution. t=%.4f' % t)
    plt.hist(v_x,bins,histtype='bar',align='mid',orientation='vertical')
    plt.ylabel('N particles')
    #plt.axis([0 26 -3 3])
    plt.subplot(3,1,2)
    plt.hist(v_y,bins,histtype='bar',align='mid',orientation='vertical')
    plt.ylabel('N particles')
    plt.subplot(3,1,3)
    plt.hist(v_z,bins,histtype='bar',align='mid',orientation='vertical')
    plt.ylabel('N particles')
    frame = '_tmp%03d.png' % i
    fig.savefig(frame)
    files.append(frame)
        
        
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 \
  -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.git")



print '---------------------------------------------------------------'
print 'time= %.4f ' % t
print 'number of particles= %g' % N 
