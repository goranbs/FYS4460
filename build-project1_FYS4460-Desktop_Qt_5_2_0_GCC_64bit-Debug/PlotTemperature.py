'''
PlotTemperature.py - project1_FYS4460

Should read output file from main.cpp and plot temperature as a
function of time.
'''

from scitools.std import linspace,zeros,pi
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
bins0 = []
bins1 = []
bins2 = []
bins3 = []
bins4 = []
bins5 = []
bins6 = []
bins7 = []
bins8 = []
bins9 = []
bins10 = []
bins11 = []
bins12 = []
bins13 = []
bins14 = []
bins15 = []
for line in file:
    Temp, time, Ek, Ep, Pressure, rmsq, bin0,bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,bin11,bin12,bin13,bin14,bin15 = line.split() # split on whitespace
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
    bins0.append(int(bin0))
    bins1.append(int(bin1))
    bins2.append(int(bin2))
    bins3.append(int(bin3))
    bins4.append(int(bin4))
    bins5.append(int(bin5))
    bins6.append(int(bin6))
    bins7.append(int(bin7))
    bins8.append(int(bin8))
    bins9.append(int(bin9))
    bins10.append(int(bin10))
    bins11.append(int(bin11))
    bins12.append(int(bin12))
    bins13.append(int(bin13))
    bins14.append(int(bin14))
    bins15.append(int(bin15))

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


binz = zeros(16)
for i in range(len(bin0)):
    binz[0] += bins0[i]
    binz[1] += bins1[i]
    binz[2] += bins2[i]
    binz[3] += bins3[i]
    binz[4] += bins4[i]
    binz[5] += bins5[i]
    binz[6] += bins6[i]
    binz[7] += bins7[i]
    binz[8] += bins8[i]
    binz[9] += bins9[i]
    binz[10] += bins10[i]
    binz[11] += bins11[i]
    binz[12] += bins12[i]
    binz[13] += bins13[i]
    binz[14] += bins14[i]
    binz[15] += bins15[i]

binz[:] = binz[:]/len(bin0)

radius = linspace(0.5,2.0,16)
volumeR = zeros(16)
volumeR[0] = (4./3)*pi*0.5**3
for ii in range(16-1):
    i = ii+1
    volumeR[i] = (4./3)*pi*(radius[i]**3 - radius[ii]**3)

rho = 1.0856

for i in range(16):
    binz[i] = binz[i]/(rho*volumeR[i])

plt.figure()
plt.plot(radius,binz,'r-*')

plt.show(True)
