#!/usr/bin/env python
# Creates an animation of the ramsey sequence with circular RF

# Some initial parameters
W_VAL = 20   #[rad s^-1]
W0_VAL = 20  #[rad s^-1]    Static field strength
WC_VAL = 3.14159265359/2 #[rad s^-1]   Rotating field strength

# Time parameters
PULSE_1_TIME = 1        # [seconds]
PULSE_2_TIME = 1        # [seconds]
RK_STEP = 0.001       # [seconds]
PRECESS_TIME = 0      # [seconds]
PHI_VAL_1 = 0  #[rad]        RF pulse inital phase for first pulse

# Initial "orientation" of neutron
# Complex and real parts of the ket
A_REAL_INIT = 1
A_COMP_INIT = 0
B_REAL_INIT = 0
B_COMP_INIT = 0


def main():
    from decimal import Decimal     # For error handling
    import numpy as np

    # Initialize stuff
    n = 4     # number of equations in the vector
    dt = RK_STEP
    t1Range = np.arange(0, PULSE_1_TIME + dt, dt)
    tPrecessRange = np.arange(0, PRECESS_TIME, dt)
    t2Range =  np.arange(0, PULSE_2_TIME + dt, dt)
    u0 = np.zeros(n)
    u0[0] = A_REAL_INIT
    u0[1] = A_COMP_INIT
    u0[2] = B_REAL_INIT
    u0[3] = B_COMP_INIT

    zProb = []
    xProb = []
    yProb = []

    # Error handling
    if ( (Decimal(str(PULSE_1_TIME)) % Decimal(str(RK_STEP)) != 0) \
        or (Decimal(str(PULSE_2_TIME)) % Decimal(str(RK_STEP)) != 0)):
        print("Error: Pulse time resolution too small for RK integrator")
        return

    # First pulse
    for t0 in t1Range:
        # odds of measuring spin along z, x, and y
        # derived from eqs 3.22 - 3.24 in May's nEDM thesis
        # u0[0] = Re[c], u0[1] = Im[c], u0[2] = Re[d], u0[3] = Im[d]

        zProb.append(u0[0]*u0[0] + u0[1]*u0[1])
        xProb.append(1/2 + u0[0]*u0[2] + u0[1]*u0[3])
        yProb.append(1/2 + u0[1]*u0[2] - u0[3]*u0[0])

        # take one RK step
        u0 = rkStep ( t0, n, u0, dt, PHI_VAL_1, spinor)

    # Larmor precess
    if len(tPrecessRange) >= 3000 :
        # Append some values to x, y, zProb to make the animation look nice
        # Pick the first and last 1000 values, let's say
        for t0 in tPrecessRange[:2000]:
            uTemp = larmor(u0, t0, W0_VAL, n)
            zProb.append(uTemp[0]*uTemp[0] + uTemp[1]*uTemp[1])
            xProb.append(1/2 + uTemp[0]*uTemp[2] + uTemp[1]*uTemp[3])
            yProb.append(1/2 + uTemp[1]*uTemp[2] - uTemp[3]*uTemp[0])

        for t0 in tPrecessRange[-1000:]:
            uTemp = larmor(u0, t0, W0_VAL, n)
            zProb.append(uTemp[0]*uTemp[0] + uTemp[1]*uTemp[1])
            xProb.append(1/2 + uTemp[0]*uTemp[2] + uTemp[1]*uTemp[3])
            yProb.append(1/2 + uTemp[1]*uTemp[2] - uTemp[3]*uTemp[0])


    # Actual larmor precession
    u0 = larmor(u0, PRECESS_TIME, W0_VAL, n)

    # Second pulse for t0 in t1Range:
    phiVal2 = W_VAL*(PULSE_1_TIME) + PHI_VAL_1 + W_VAL*PRECESS_TIME
    for t0 in t2Range:
        zProb.append(u0[0]*u0[0] + u0[1]*u0[1])
        xProb.append(1/2 + u0[0]*u0[2] + u0[1]*u0[3])
        yProb.append(1/2 + u0[1]*u0[2] - u0[3]*u0[0])

        # take one RK step
        u0 = rkStep ( t0, n, u0, dt, phiVal2, spinor)

    print("Ending Z component: " + str(zProb[-1]))
    plotStuff(xProb, yProb, zProb)

    return

def plotStuff(xProb, yProb, zProb):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.animation as animation

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title('Odds of measuring spin up along x, y, z')
    # ax.scatter(xProb, yProb, zProb)
    ax.set_xlim3d(0,1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(0,1)
    ax.set_xlabel('P(x)')
    ax.set_ylabel('P(y)')
    ax.set_zlabel('P(z)')
    ax.view_init(30,220) # So that the viewing angle looks ok

    data =  np.array([xProb[::10], yProb[::10], zProb[::10]])
    line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], ".")

    # Animate and save
    print("Making animation...")
    fpsVal = 30
    nArgs = len(xProb[::10])
    ani = animation.FuncAnimation(fig, update, nArgs, fargs=(data, line), interval=1, blit=False)
    ani.save('rkCircular.mp4', fps=fpsVal)
    print("Saved to rkCircular.mp4")

    # plt.show()

    return

# no clue how this works, but it does
def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])

def spinor(t, n, u, phi):
# Modified right hand side of eq A.1 - A.4 in May's nEDM thesis
# using eq 3.28, 3.29 as a basis
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im(b)
    import numpy as np
    x = W_VAL*t + phi
    value = np.array ( [ 1/2*(W0_VAL*u[1] + WC_VAL*np.cos(x)*u[3]) - WC_VAL/2*u[2]*np.sin(x), \
                       1/2*(-W0_VAL*u[0] - WC_VAL*np.cos(x)*u[2]) - WC_VAL/2*u[3]*np.sin(x), \
                       1/2*(-W0_VAL*u[3] + WC_VAL*np.cos(x)*u[1]) + WC_VAL/2*u[0]*np.sin(x), \
                       1/2*(W0_VAL*u[2] - WC_VAL*np.cos(x)*u[0]) + WC_VAL/2*u[1]*np.sin(x) ])

    return value

def larmor(u0, dt, w0, n):
# Analytical larmor precession period for some time dt
# Takes in an initial ket u0=(Re(a),Im(a),Re(b),Im(b))
# and returns a final ket ketFinal. Requires n = 4 elements in u0
    import numpy as np
    if ( n != 4 ):
        print("Please set n = 4")
        return

    x = dt*w0/2
    ketFinal = np.array ( [u0[0]*np.cos(x) + u0[1]*np.sin(x), \
                        u0[1]*np.cos(x) - u0[0]*np.sin(x), \
                        u0[2]*np.cos(x) - u0[3]*np.sin(x), \
                        u0[3]*np.cos(x) + u0[2]*np.sin(x)  ])
    return ketFinal

def rkStep ( t0, m, u0, dt, phi, f):
# Takes one runge kutta step
# adapted from rk4 library to accomodate an f with frequency and phase
    import numpy as np
#
#  Get four sample values of the derivative.
#
    f0 = f ( t0, m, u0, phi )

    t1 = t0 + dt / 2.0
    u1 = np.zeros ( m )
    u1[0:m] = u0[0:m] + dt * f0[0:m] / 2.0
    f1 = f ( t1, m, u1, phi )

    t2 = t0 + dt / 2.0
    u2 = np.zeros ( m )
    u2[0:m] = u0[0:m] + dt * f1[0:m] / 2.0
    f2 = f ( t2, m, u2, phi )

    t3 = t0 + dt
    u3 = np.zeros ( m )
    u3[0:m] = u0[0:m] + dt * f2[0:m]
    f3 = f ( t3, m, u3, phi )
    #
    #  Combine them to estimate the solution U at time T1.
    #
    u = np.zeros ( m )
    u[0:m] = u0[0:m] + ( dt / 6.0 ) * ( \
            f0[0:m] \
    + 2.0 * f1[0:m] \
    + 2.0 * f2[0:m] \
    +       f3[0:m] )

    return u

if ( __name__ == '__main__' ):
    main()
