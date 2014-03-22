'''
PlotTemperature.py - project1_FYS4460

Should read output file from main.cpp and plot temperature as a
function of time.
'''

from scitools.std import linspace,zeros,pi
filename = 'temperatures.txt'
file = open(filename,'r')

nbins = int(file.readline()) # first line is number of bins
rho = float(file.readline()) # second is density of system
N = int(file.readline())     # third line is number of particles in system
file.readline() # read fourth and fifth line
file.readline()

E_kin = []   # kinetic energy
E_pot = []   # potetial energy
E_tot = []   # total energy of the system 
T = []       # Temperataure
P = []       # Pressure
r_msq_t = [] # Mean square displacement
t = []       # time
nsy = []     # not sure yet :-)

binz = zeros(nbins)
timeiterations = 0
for line in file:
    values = line.split()

    Temp = values[0]
    time = values[1]
    Ek = values[2]
    Ep = values[3]
    Pressure = values[4]
    rmsq = values[5]
    allbins = values[6:]

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

    timeiterations += 1
    for i in range(len(allbins)):
        binz[i] += (int(allbins[i]))
        


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


#######################################3
print timeiterations
l_binz = len(binz)
print l_binz

radius = linspace(0,(l_binz+1)*0.1,l_binz)

#plt.figure()
#plt.plot(radius,binz)


volumeR = zeros(l_binz)

for i in range(l_binz-1):
    volumeR[i] = (4./3)*pi*(radius[i+1]**3 - radius[i]**3) # volume of shell in range [r,r+dr]

binz[1:] = binz[1:]/(timeiterations*N*rho*volumeR[i]) # mean number of particles in volume

plt.figure()
plt.plot(radius,binz,'r-*')
plt.title('Average number of particles in distance r from atom')
plt.xlabel('r [MD-units]')
plt.ylabel('Number of particles')
plt.legend(['N(r)'])

# g(r) function: probability of finding an atom in distance r from another atom.
integral = 0
for i in range(l_binz-1):
    integral += (binz[i+1] + binz[i])

integral = integral*0.5*0.1
#print integral
plt.figure()
plt.plot(radius,binz/integral,'b-*')
plt.title('Particle density in distance r from atom')
plt.xlabel('r [MD-units]')
plt.ylabel('particle density')
plt.legend(['g(r)'])


plt.show(True)


