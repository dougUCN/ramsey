#!/usr/bin/env python
#
#   4th order runge kutta on linear on nEDM in
#   a Magnetic Field with Uniform Static and Transverse
#   Circularly Oscillating Components
#
#   Douglas Wong 2/20/17

# Some initial parameters
W_VAL = 23   #[rad s^-1]
W0_VAL = 20  #[rad s^-1]    Static field strength
WC_VAL = 1.57 #[rad s^-1]   Rotating field strength
PHI_VAL = 21 #[rad]          RF pulse inital phase

# Initial "orientation" of neutron
# Complex and real parts of the ket
A_REAL_INIT = -0.82573592
A_COMP_INIT = 0.32626352
B_REAL_INIT = -0.40023091
B_COMP_INIT = -0.22699674

# Step taken by integrator and total period
MAX_TIME = 1        # [seconds]
TIME_STEP = 0.001    # [seconds]


def main():
    from rk4 import rk4vec
    from rk4 import rk4vec_test_f
    import numpy as np

    # Initialize stuff
    n = 4     # number of equations in the vector
    PI_CONST = 3.14159265359
    dt = TIME_STEP
    tRange = np.arange(0, MAX_TIME + dt, dt)
    u0 = np.zeros(n)
    u0[0] = A_REAL_INIT
    u0[1] = A_COMP_INIT
    u0[2] = B_REAL_INIT
    u0[3] = B_COMP_INIT

    zProb = []
    xProb = []
    yProb = []

    # x - y components of circularly rotating magnetic B field
    xB = []
    yB = []
    angleDiff = [] #Between semiclassical neutron orientation in x-y and the field


    for t0 in tRange:
        # odds of measuring spin along z, x, and y
        # derived from eqs 3.22 - 3.24 in May's nEDM thesis
        # u0[0] = Re[c], u0[1] = Im[c], u0[2] = Re[d], u0[3] = Im[d]

        zProb.append(u0[0]*u0[0] + u0[1]*u0[1])
        xProb.append(1/2 + u0[0]*u0[2] + u0[1]*u0[3])
        yProb.append(1/2 + u0[1]*u0[2] - u0[3]*u0[0])

        xB.append(np.cos(W_VAL*(-1*t0) + PHI_VAL))
        yB.append(np.sin(W_VAL*(-1*t0) + PHI_VAL))

        xTemp = xProb[-1] - 0.5
        yTemp = yProb[-1] - 0.5

        # Dot product to find angle between RF field and semiclassical
        # neutron orientation
        angleDiff.append(180/np.pi * np.arccos((xTemp*xB[-1] + yTemp*yB[-1]) \
                            / (np.sqrt(xTemp*xTemp + yTemp*yTemp)) \
                            / (np.sqrt(xB[-1]*xB[-1] + yB[-1]*yB[-1])) ) )

        # take one RK step
        u0 = rk4vec ( t0, n, u0, dt, spinor)
    # END FOR
    print(u0)

    plotStuff(xProb, yProb, zProb, tRange)

    return

def plotStuff(xProb, yProb, zProb, time):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig2 = plt.figure(2)
    ax = fig2.add_subplot(111, projection='3d')
    ax.plot(xProb, yProb, zProb, color="#005F73")
    ax.set_xlim3d(0,1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(0,1)
    ax.set_xlabel('P(x)')
    ax.set_ylabel('P(y)')
    ax.set_zlabel('P(z)')
    ax.view_init(12,220) # So that the viewing angle looks ok

    plt.figure(figsize=[4.8,6])
    ax1 = plt.subplot(311)
    ax1.plot(time, xProb, color="#AE2012")
    ax1.set(xlim=[0, time[-1]], ylim=[0,1], ylabel='P(x)')
    ax1.xaxis.set_tick_params(labelbottom=False)

    ax2 = plt.subplot(312, sharex=ax1, sharey=ax1)
    ax2.plot(time, yProb, color="#EE9B00")
    ax2.set(ylabel='P(y)')
    ax2.xaxis.set_tick_params(labelbottom=False)
    
    ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
    ax3.plot(time, zProb, color="#0A9396")
    ax3.set(ylabel='P(z)', xlabel='time [s]')

    plt.show()

    return

def spinor(t, n, u):
# Modified right hand side of eq A.1 - A.4 in May's nEDM thesis
# using eq 3.28, 3.29 as a basis
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im(b)
    import numpy as np
    x = W_VAL*t + PHI_VAL
    value = np.array ( [ 1/2*(W0_VAL*u[1] + WC_VAL*np.cos(x)*u[3]) - WC_VAL/2*u[2]*np.sin(x), \
                       1/2*(-W0_VAL*u[0] - WC_VAL*np.cos(x)*u[2]) - WC_VAL/2*u[3]*np.sin(x), \
                       1/2*(-W0_VAL*u[3] + WC_VAL*np.cos(x)*u[1]) + WC_VAL/2*u[0]*np.sin(x), \
                       1/2*(W0_VAL*u[2] - WC_VAL*np.cos(x)*u[0]) + WC_VAL/2*u[1]*np.sin(x) ])

    return value

if ( __name__ == '__main__' ):
    main()
