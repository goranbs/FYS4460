'''
This program should read velocity data from stateXXXX.txt file generated
by the MD program main.cpp and plot histograms of the velocity distribution
of the particles.
The particles should eventually after some time have a velocity distribution
similar to a Maxwell-Bolzmann distribution.
'''

import matplotlib.pyplot as plt
import numpy as np
import os, sys


def directory(path,extension):
  list_dir = []
  list_dir = os.listdir(path)
  count = 0
  for file in list_dir:
    if file.endswith(extension): # eg: '.txt'
      count += 1
  return count

num_files = directory('/home/goranbs/goran/FYS4460 - Percolation/build-project1_FYS4460-Desktop_Qt_5_2_0_GCC_64bit-Debug/','.txt')

print 'Directory PATH:/home/goranbs/goran/FYS4460 - Percolation/build-project1_FYS4460-Desktop_Qt_5_2_0_GCC_64bit-Debug/'
print 'Number of files in directory= %g' % num_files

filename = 'state%04g.txt' % 0      # filename = 'state0000.txt' initial file
file = open(filename,'r')

N = int(file.readline())            # number of atoms
secondline = file.readline()
somestring,t = secondline.split()
t = float(t)                        # time in simulation

files = []
num_files_used = 0;
for i in range(int((num_files-2)/4.0)):
  
  i = i*4 # every fourth file
  num_files_used += 1

  filename = 'state%04g.txt' % i      # filename = 'state000.txt'
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
    v_x.append(abs(float(line[4])))
    v_y.append(abs(float(line[5])))
    v_z.append(abs(float(line[6])))
    j += 1
        
  if j == N_atoms:
    print 'reading file nr %g OK. Processing image...' % i
  else:
    print ' j != N, j=%g, N=%g' % j,N
    
  # plotting    
  
  bins = 10
  fontsize = 'medium'
  fig = plt.figure(figsize=(5,5))
  #plt.setsize(fontsize)
  plt.subplot(3,1,1)
  plt.title('velocity distribution. t=%.2f [fs]' % t)
  plt.hist(v_x,bins,histtype='bar',align='mid',orientation='vertical')
  plt.ylabel('N particles')
  #plt.axis([0 26 -3 3])
  plt.subplot(3,1,2)
  plt.hist(v_y,bins,histtype='bar',align='mid',orientation='vertical')
  plt.ylabel('N particles')
  plt.subplot(3,1,3)
  plt.hist(v_z,bins,histtype='bar',align='mid',orientation='vertical')
  plt.ylabel('N particles')
  frame = '_tmp%04d.png' % num_files_used
  fig.savefig(frame)
  files.append(frame)
        
        
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 \
  -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.git")



print '---------------------------------------------------------------'
print 'time= %.4f ' % t
print 'number of particles= %g' % N 
print 'number of files counted= %g' % num_files-1
