# -*- coding: utf-8 -*-
"""
Created on Thu May  1 13:48:45 2014

@author: goranbs
"""

# i) Percolation project 3 - FYS4460

# Read inputfile PI_lattices.dat, find p for PI = 0.3 and PI = 0.8 by interpolation
# Plot p as function of latticesizes L.


from scitools.std import *
from numpy import interp, zeros, polyfit, polyval

filename = 'PI_lattices.dat' # found from the matlab program i_find_pPI.m 
file = open(filename,'r') # read access

firstline = file.readline()
secondline = file.readline()

L = [float(x) for x in secondline.split()]
print L

len_L = len(L)

PI_1 = []
PI_2 = []
PI_3 = []
PI_4 = []
PI_5 = []
PI_6 = []
p = []
for line in file:
    pi1,pi2,pi3,pi4,pi5,pi6,prob = [float(x) for x in line.split()]
    PI_1.append(pi1)
    PI_2.append(pi2)
    PI_3.append(pi3)
    PI_4.append(pi4)
    PI_5.append(pi5)
    PI_6.append(pi6)
    p.append(prob)

len_p = len(p)

PI = zeros((len_p,len_L))

for i in range(len_p):
    PI[i,0] = PI_1[i]
    PI[i,1] = PI_2[i]
    PI[i,2] = PI_3[i]
    PI[i,3] = PI_4[i]
    PI[i,4] = PI_5[i]
    PI[i,5] = PI_6[i]
    
# testplotting:

plot(p,PI_1,'r-')    
hold('on')
plot(p,PI_2,'b-')
plot(p,PI_3,'y-')
plot(p,PI_4,'g-')
plot(p,PI_5,'m-')
plot(p,PI_6,'bk-')
legend('L1','L2','L3','L4','L5','L6')
title('PI(p)')
hold('off')

# find p for PI = 0.3 and PI = 0.8


x1 = 0.3
x2 = 0.8

p_values_x1 = []
p_values_x2 = []

for i in range(len_L):
    a = interp(x1,PI[:,i],p)
    b = interp(x2,PI[:,i],p)
    
    p_values_x1.append(a)
    p_values_x2.append(b)
    
    
    
# Plotting
import matplotlib.pyplot as plt
Fontsize = 20
name = 'i_find_p_for_PI.png'
name2 = 'j_estimate_nu.png'
name3 = 'k_estimate_pc.png'

fig = plt.figure()
plt.plot(L,p_values_x1,'b-*')
plt.hold(True)
plt.plot(L,p_values_x2,'r-*')
plt.title(r'$p_{\Pi =x}$ as function of lattice size',fontsize=Fontsize)
plt.xlabel('Lattice sizes',fontsize=Fontsize)
plt.ylabel('cutoffprobability',fontsize=Fontsize)
plt.legend((r'$p_{\Pi =0.3}$',r'$p_{\Pi =0.8}$'),fontsize=Fontsize)
plt.savefig(name)
plt.hold(False)


log_L = []
log_pp = []
for i in range(len_L):
    x = log10(p_values_x2[i] - p_values_x1[i])
    log_pp.append(x)
    log_L.append(log10(L[i]))

pp = polyfit(log_L,log_pp,1)
f = polyval(pp,log_L)



#print pp[0], pp[1]
nu = -1/pp[0]
error = [(log_pp[i] - f[i])**2 for i in range(len(log_pp))]
stdev = mean(error)

fig2 = plt.figure()
plt.plot(log_L,log_pp,'b-*')
plt.hold(True)
plt.plot(log_L,f,'r--')
plt.title(r'Estimate on $ \nu $ = %.4f $ \pm $ %.5f' % (nu,stdev),fontsize=Fontsize)
plt.xlabel(r'$log(L)$', fontsize=Fontsize)
plt.ylabel(r'$log(p_{\Pi =0.8} - p_{\Pi =0.3})$', fontsize=Fontsize)
plt.savefig(name2)
plt.hold(False)

nu_exact = 0.8751
L_exp = [L[i]**(-1/nu_exact) for i in range(len(L))]

fig3 = plt.figure()
plt.plot(L_exp, p_values_x1, 'b-o')
plt.hold(True)
plt.plot(L_exp, p_values_x2, 'r-o')
plt.title(r'Estimate $ p_c $', fontsize=Fontsize)
plt.xlabel(r'$ L^{-1/\nu} $', fontsize=Fontsize)
plt.ylabel(r'$ p_{\Pi = x} $', fontsize=Fontsize)
plt.legend((r'$ p_{\Pi = 0.3} $',r'$ p_{\Pi = 0.8} $'), loc='lower left', fontsize=Fontsize)
plt.savefig(name3)
plt.hold(False)

plt.show(True)