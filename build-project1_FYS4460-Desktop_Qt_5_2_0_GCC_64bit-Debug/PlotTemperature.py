'''
PlotTemperature.py - project1_FYS4460

Should read output file from main.cpp and plot temperature as a
function of time.
'''


filename = 'temperatures.txt'
file = open(filename,'r')

file.readline() # read fist and second line
file.readline()
E_kin = [] # kinetic energy
T = [] # Temperataure
t = [] # time
for line in file:
    Temp, time, Ek = line.split() # split on whitespace
    T.append(float(Temp))
    t.append(float(time))
    E_kin.append(float(Ek))

import matplotlib.pyplot as plt

plt.figure()
plt.plot(t,T)
plt.title('Temperature of system as funciton of time')
plt.xlabel('time [fs]')
plt.ylabel('Temperature [K]')
plt.legend('T(t)')

plt.figure()
plt.plot(t,E_kin)
plt.title('Kinetic energy as a function of time')
plt.xlabel('time [fs]')
plt.ylabel('Energy [eV]')
plt.legend('Ek(t)')
plt.show(True)
