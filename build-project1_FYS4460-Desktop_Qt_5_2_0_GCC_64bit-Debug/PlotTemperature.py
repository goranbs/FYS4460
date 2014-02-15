'''
PlotTemperature.py - project1_FYS4460

Should read output file from main.cpp and plot temperature as a
function of time.
'''


filename = 'temperatures.txt'
file = open(filename,'r')

file.readline() # read fist and second line
file.readline()
T = [] # Temperataure
t = [] # time
for line in file:
    Temp, time = line.split() # split on whitespace
    T.append(float(Temp))
    t.append(float(time))

import matplotlib.pyplot as plt

plt.figure()
plt.plot(t,T)
plt.title('Temperature of system as funciton of time')
plt.xlabel('unitless time')
plt.ylabel('unitless temperature')
plt.legend('T(t)')
plt.show(True)
