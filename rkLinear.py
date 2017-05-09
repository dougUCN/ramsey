#!/usr/bin/env python
#
#   4th order runge kutta on linear on nEDM in
#   a Magnetic Field with Uniform Static and Transverse
#   Linearly Oscillating Components
#
#   Douglas Wong 2/8/17

# Some initial parameters
W_VAL = 20   #[rad s^-1]
W0_VAL = 20  #[rad s^-1]
WL_VAL = 3.1415 #[rad s^-1]
PHI_VAL = 0  #[rad]

# Initial "orientation" of neutron
A_INIT = 1
B_INIT = 0

# Step taken by integrator and total period
MAX_TIME = 2        # [seconds]
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

    plt.figure(1)
    plt.plot(time, zProb)
    plt.title('Odds of measuring spin up along z')
    plt.xlabel('time [s]')
    plt.ylabel('P(z)')
    plt.axis([0,2*np.pi,0,1])

    fig2 = plt.figure(2)
    ax = fig2.add_subplot(111, projection='3d')
    plt.title('Odds of measuring spin along x, y, z')
    ax.scatter(xProb, yProb, zProb)
    ax.set_xlabel('P(x)')
    ax.set_ylabel('P(y)')
    ax.set_zlabel('P(z)')
    ax.view_init(30,220) # So that the viewing angle looks ok

    plt.figure(3)
    plt.plot(time, xProb)
    plt.title('Odds of measuring spin up along x')
    plt.xlabel('time [s]')
    plt.ylabel('P(x)')
    plt.axis([0,2*np.pi,0,1])

    plt.figure(4)
    plt.plot(time, yProb)
    plt.title('Odds of measuring spin up along y')
    plt.xlabel('time [s]')
    plt.ylabel('P(y)')
    plt.axis([0,2*np.pi,0,1])

    plt.show()

    return

if ( __name__ == '__main__' ):
    main()
