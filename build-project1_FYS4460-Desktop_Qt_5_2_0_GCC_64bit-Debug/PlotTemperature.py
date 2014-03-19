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


###################################################################
# Mean Temperature and mean pressure

from scitools.std import sqrt

def Temperature_and_Pressure():
    meanT = 0
    meanP = 0
    mean_squareT = 0
    mean_squareP = 0
    # begin to sample after 200 timesteps:
    N = len(T) -200
    for j in range(N):
        i = j+200
        meanT += T[i]
        meanP += P[i]
        mean_squareT += T[i]*T[i]
        mean_squareP += P[i]*P[i]
        
    meanT = meanT/N
    meanP = meanP/N

    sT = sqrt(mean_squareT/N - meanT*meanT)
    sP = sqrt(mean_squareP/N - meanP*meanP)

    return meanT,meanP,sT,sP


[meanT,meanP,sT,sP] = Temperature_and_Pressure()

print 'Mean temp  = %.2f pm %.3f   MD-units ' % (meanT,sT)
print 'Mean press = %.2f pm %.3f   MD-units ' % (meanP,sP)
 

###################################################################
# Plotting

D = r_msq_t[-1]/(t[-1]*6)

# print 'Estimated diffusion constant D=%.1f ' % (D)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(t,T)
plt.title('Temperature of system as funciton of time. T= %.2f pm %.3f   MD-units' % (meanT,sT))
plt.xlabel('time [MD]')
plt.ylabel('Temperature [MD]')
plt.legend('T(t)')

plt.figure()
plt.plot(t,E_kin, 'r-')
plt.hold(True)
plt.plot(t,E_pot, 'm-')
plt.plot(t,E_tot, 'b-*')
plt.title('Energy as a function of time')
plt.xlabel('time [MD]')
plt.ylabel('Energy [MD]')
plt.legend(('Ek(t)','U(t)', 'E_{tot}(t)'), loc='lower right')
plt.hold(False)

plt.figure()
plt.plot(t,P,'r-')
plt.title('pressure of system as function of time. P= %.1f pm %.3f  MD-units' % (meanP,sP))
plt.xlabel('time [MD]')
plt.ylabel('pressure [MD]')
plt.legend(('P(t)'), loc='lower right')

plt.figure()
plt.plot(t,r_msq_t,'b-d')
plt.title('mean square displacement of particles in system')
plt.ylabel('time [MD]')
#plt.xlabel('time [fs]')
#plt.ylabel('msq [m^2]')
plt.ylabel('msq [MD]')
plt.legend(('MSQ(t)'), loc='lower right')

plt.figure()
plt.plot(t,nsy,'r-')
plt.title('Diffusion  D = %.1f ' % D)
plt.xlabel('time [MD]')
#plt.xlabel('time [fs]')
plt.ylabel('msq/t [MD]')
#plt.ylabel('msq/t [m^2/t]')
plt.legend(('Diffusion constant'), loc='lower right')




plt.show(True)
