'''
PlotTemperature.py - project1_FYS4460

Should read output file from main.cpp and plot temperature as a
function of time.
'''


filename = 'temperatures.txt'
file = open(filename,'r')

file.readline() # read fist and second line
file.readline()
E_kin = []   # kinetic energy
E_pot = []   # potetial energy
E_tot = []   # total energy of the system 
T = []       # Temperataure
P = []       # Pressure
r_msq_t = [] # Mean square displacement
t = []       # time
nsy = []     # not sure yet :-)
for line in file:
    Temp, time, Ek, Ep, Pressure, rmsq = line.split() # split on whitespace
    T.append(float(Temp))
    t.append(float(time))
    E_kin.append(float(Ek))
    E_pot.append(float(Ep))
    E_tot.append(float(Ek) + float(Ep))
    P.append(float(Pressure))
    r_msq_t.append(float(rmsq))
    if (float(time) == 0):
        nsy.append(0)
    else:
        nsy.append(float(rmsq)/(float(time)*6))

import matplotlib.pyplot as plt

plt.figure()
plt.plot(t,T)
plt.title('Temperature of system as funciton of time')
plt.xlabel('time [fs]')
plt.ylabel('Temperature [K]')
plt.legend('T(t)')

plt.figure()
plt.plot(t,E_kin, 'r-')
plt.hold(True)
plt.plot(t,E_pot, 'm-')
plt.plot(t,E_tot, 'b-*')
plt.title('Energy as a function of time')
plt.xlabel('time [fs]')
plt.ylabel('Energy [eV]')
plt.legend(('Ek(t)','U(t)', 'E_{tot}(t)'), loc='lower right')
plt.hold(False)

plt.figure()
plt.plot(t,P,'r-')
plt.title('pressure of system as function of time')
plt.xlabel('time [fs]')
plt.ylabel('pressure [N/m^2]')
plt.legend(('P(t)'), loc='lower right')

plt.figure()
plt.plot(t,r_msq_t,'b-d')
plt.title('mean square displacement of particles in system')
plt.xlabel('time [fs]')
plt.ylabel('msq [m^2]')
plt.legend(('MSQ(t)'), loc='lower right')

plt.figure()
plt.plot(t,nsy,'r-')
plt.title('Diffusion constant of particles in system')
plt.xlabel('time [fs]')
plt.ylabel('msq/t [m^2/t]')
plt.legend(('6Dt'), loc='lower right')

D = r_msq_t[-1]/(t[-1]*6)

print 'Estimated diffusion constant D=%.1f ' % (D)

#plt.figure()
#plt.plot(t,E_kin, 'r-')
#plt.title('Kinetic energy as a function of time')
#plt.xlabel('time [fs]')
#lt.ylabel('Energy [eV]')
#lt.legend('Ek(t)', loc='upper right')

#plt.figure()
#plt.plot(t,E_pot, 'r-')
#plt.title('Potential energy as a function of time')
#plt.xlabel('time [fs]')
#plt.ylabel('Energy [eV]')
#plt.legend('U(t)', loc='upper right')

#plt.figure()
#plt.plot(t,E_tot, 'r-')
#plt.title('Total energy as a function of time')
#plt.xlabel('time [fs]')
#plt.ylabel('Energy [eV]')
#lt.legend('E(t)', loc='upper right')

plt.show(True)
