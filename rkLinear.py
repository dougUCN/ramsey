#!/usr/bin/env python
#
#   4th order runge kutta on linear on nEDM in
#   a Magnetic Field with Uniform Static and Transverse
#   Linearly Oscillating Components
#
#   Douglas Wong 2/8/17

# Some initial parameters
W_VAL = 20    #[rad s^-1]
W0_VAL = 20  #[rad s^-1]
WL_VAL = 1.57/2 #[rad s^-1]
PHI_VAL = 0  #[rad]

# Initial "orientation" of neutron
A_INIT = 1
B_INIT = 0

# Step taken by integrator and total period
MAX_TIME = 8        # [seconds]
TIME_STEP = 0.001    # [seconds]

def main():
    from rk4 import rk4vec
    from rk4 import rk4vec_test_f
    import numpy as np
    import matplotlib.pyplot as plt

    # Initialize stuff
    n = 4     # number of equations in the vector
    dt = TIME_STEP
    tRange = np.arange(0, MAX_TIME + dt, dt)

    u0 = np.zeros ( n )
    u0[0] = A_INIT
    u0[2] = B_INIT

    zProb = []
    xProb = []
    yProb = []

    for t0 in tRange:
        # odds of measuring spin along z, x, and y
        # u0[0] = Re[a], u0[1] = Im[a], u0[2] = Re[b], u0[3] = Im[b]

        zProb.append(u0[0]*u0[0] + u0[1]*u0[1])
        xProb.append(1/2 + u0[0]*u0[2] + u0[1]*u0[3])
        yProb.append(1/2 + u0[1]*u0[2] - u0[3]*u0[0])

        # Takes one RK step
        u0 = rk4vec ( t0, n, u0, dt, spinor)
    #  END FOR

    plotStuff(xProb, yProb, zProb, tRange)

    return

def spinor(t, n, u):
# Right hand side of eq A.1 - A.4 in May nEDM thesis
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im b
    import numpy as np
    w = W_VAL   #[rad s^-1]
    w0 = W0_VAL  #[rad s^-1]
    wl = WL_VAL #[rad s^-1]
    phi = PHI_VAL   #[rad]
    value = np.array ( [ 1/2*(w0*u[1] + wl*np.cos(w*t + phi)*u[3]), \
                       1/2*(-w0*u[0] - wl*np.cos(w*t + phi)*u[2]), \
                       1/2*(-w0*u[3] + wl*np.cos(w*t + phi)*u[1]), \
                       1/2*(w0*u[2] - wl*np.cos(w*t + phi)*u[0])])

    return value

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
    ax.view_init(30,220) # So that the viewing angle looks ok

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
    ins = ax3.inset_axes([0.7,0.5, 0.25, 0.4]) # relative position on x, y, and relative x y box height?
    ins.plot(time, zProb, color="#0A9396")
    ins.set_xlim(4, 4.1)
    ins.set_ylim(0.45, 0.55)
    ins.indicate_inset_zoom(ins, edgecolor="black")

    plt.show()

    return

if ( __name__ == '__main__' ):
    main()
