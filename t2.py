#!/usr/bin/env python
#Estimates t2 time from fixed neutrons in a gradient
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit

def fit(t, c, b, w, phi):
    # Function for fitting <Sx>
    return np.exp(-(t - c) / b ) * np.sin(w*t + phi)

# Guess for parameters to fit
guess = [30, 50,183, 807]

def larmor(u0, dt, w0):
# Analytical larmor precession period for some time dt
# Takes in an initial ket u0=(Re(a),Im(a),Re(b),Im(b))
# and returns a final ket ketFinal
    x = dt*w0/2
    ketFinal = np.array ( [u0[0]*np.cos(x) + u0[1]*np.sin(x), \
                        u0[1]*np.cos(x) - u0[0]*np.sin(x), \
                        u0[2]*np.cos(x) - u0[3]*np.sin(x), \
                        u0[3]*np.cos(x) + u0[2]*np.sin(x)  ])
    return ketFinal

def getXProb(u):
    return 0.5 + u[0]*u[2] + u[1]*u[3]

def getYProb(u):
    return 0.5 + u[1]*u[2] - u[3]*u[0]

def getZProb(u):
    return u[0]*u[0] + u[1]*u[1]

def main():
    t2Times = np.arange(0,30,.001)
    u0 = [1/np.sqrt(2),0,1/np.sqrt(2),0]
    b0 = 1e-6 * 1.83247172 * 10**8  # [rad/s]
    dBdZ = 1e-7 * 1.83247172 * 10**8 # [rad/(s m)]
    height = 0.1 # [m]
    wRange = np.arange(b0 - dBdZ*height/2, b0 + dBdZ*height/2, 0.1)

    sx = []
    sy = []
    for w in tqdm(wRange):
        for t2 in tqdm(t2Times):
            ket = larmor(u0, t2, w)
            sx.append( (getXProb(ket)-0.5)*2 )
            sy.append( (getYProb(ket)-0.5)*2 )

    sx = np.array(sx).reshape(len(wRange), len(t2Times))
    sy = np.array(sy).reshape(len(wRange), len(t2Times))

    # sxAv =[]
    # syAv =[]
    # for sxSet, sySet in zip(sx.transpose(),sy.transpose()):
    #     sxAv.append( np.average(sxSet) )
    #     syAv.append( np.average(sySet))
    #
    # fig = plt.figure()
    # plt.xlabel('t [sec]')
    # plt.ylabel('<Sx>')
    # plt.plot(t2Times, sxAv,color='C0')
    #
    # fig = plt.figure()
    # plt.xlabel('t [sec]')
    # plt.ylabel('<Sy>')
    # plt.plot(t2Times, syAv,color='C0')

    # Phi calculation: Make phase continuous
    phase = []
    for sxTemp, syTemp in zip(sx, sy):
        phaseTemp = np.arctan2(syTemp, sxTemp)
        discont = np.where(np.abs(np.diff(phaseTemp))> 6.0)[0] + 1  # pi -> -pi discontinuities from atan2
        discont = np.append(discont, len(phaseTemp)-1)
        counter = 1
        for i, j in zip(discont[:-1],discont[1:]):
            phaseTemp[i:j] -= 2*counter*np.pi           # MIGHT NEED TO BE += or -=
            counter += 1
        phase.append( phaseTemp )

    phaseStd = []
    phaseAv = []
    for phaseSet in np.transpose(phase):
        phaseStd.append( np.std(phaseSet) )
        phaseAv.append( np.average(phaseSet) )

    fig = plt.figure()
    for neutron in phase:
        plt.plot(t2Times[:-1], neutron[:-1],color='C0')
    plt.grid(True)
    plt.xlabel('t [s]')
    plt.ylabel('Phase [rad]')

    fig = plt.figure()
    plt.plot(t2Times[:-1], phaseStd[:-1])
    plt.grid(True)
    plt.xlabel('t')
    plt.ylabel('$\sigma$')

    fig = plt.figure()
    plt.plot(t2Times, np.array(phaseStd)/np.array(phaseAv))
    plt.grid(True)
    plt.xlabel('t')
    plt.ylabel('$\sigma$/<$\phi$>')

    plt.show()



    return

if ( __name__ == '__main__' ):
    main()
