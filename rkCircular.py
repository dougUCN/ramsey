#!/usr/bin/env python
#
#   4th order runge kutta on linear on nEDM in
#   a Magnetic Field with Uniform Static and Transverse
#   Circularly Oscillating Components
#
#   Douglas Wong 2/20/17

W_DRIVE = 20   #[rad s^-1]
W0 = 20  #[rad s^-1]
WC = 1.57 #[rad s^-1]
PHI = 3.14159/2  #[rad]

# Initial orientation of neutron
A_INIT = 0.7071
B_INIT = 0.7071

def main():
    from rk4 import rk4vec
    from rk4 import rk4vec_test_f
    import numpy as np

    # Initialize stuff
    n = 4     # number of equations in the vector
    dt = 0.001
    tmax = 1
    t0 = 0

    u0 = np.zeros(n)
    u0[0] = A_INIT
    u0[1] = 0
    u0[2] = B_INIT
    u0[3] = 0

    time = [0]
    zProb = []
    xProb = []
    yProb = []

    i = 0
    while (True):
        # print ( '  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g' \
        #         % ( i, time[i], u0[0], u0[1], u0[2], u0[3]))

        # odds of measuring spin along z, x, and y
        # derived from eqs 3.22 - 3.24 in May's nEDM thesis
        # u0[0] = Re[c], u0[1] = Im[c], u0[2] = Re[d], u0[3] = Im[d]

        # USE THIS ONLY FOR SPINORTEST
        # These come from explicity using eq 3.30-3.31 in the May thesis
        # temp = (W_DRIVE*t0 + PHI)
        #
        # zProb.append(np.power(u0[0], 2) + np.power(u0[1], 2))
        #
        # xProb.append(1/2 + (u0[0]*u0[2] + u0[1]*u0[3])*np.cos(temp)\
        #                 + (u0[2]*u0[1] - u0[0]*u0[3])*np.sin(temp))
        # yProb.append(1/2 + (u0[2]*u0[1] - u0[0]*u0[3])*np.cos(temp)\
        #                 - (u0[0]*u0[2] + u0[1]*u0[3])*np.sin(temp))

        # USE THIS FOR SPINOR
        zProb.append(np.power(u0[0], 2) + np.power(u0[1], 2))
        xProb.append(1/2 + u0[0]*u0[2] + u0[1]*u0[3])
        yProb.append(1/2 + u0[1]*u0[2] - u0[3]*u0[0])

        if ( tmax <= t0 ):
            break

        i = i + 1

        # take one RK step
        t1 = t0 + dt
        u1 = rk4vec ( t0, n, u0, dt, spinor)

        # shift data
        t0 = t1
        u0 = u1.copy ( )
        time.append(t1)
  #  END WHILE

    plotStuff(xProb, yProb, zProb, time)

    print ( '' )
    print ( 'RKLINEAR:' )
    print ( '  Normal end of execution.' )
    return

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
    plt.title('Odds of measuring spin up along x, y, z')
    ax.scatter(xProb, yProb, zProb)
    ax.set_xlim3d(0,1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(0,1)
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

def spinorTest(t, n, u):
# Modified form of RHS of eq 3.32 - 3.33 in May's nEDM thesis,
# following what he did in Appendix A to linearly polarized equations
# u[0] = Re(c), u[1] = Im(c), u[2] = Re(c), u(3) = Im c
    import numpy as np
    value = np.array ( [ 1/2*(-u[1]*(W_DRIVE - W0) + u[3]*WC), \
                       1/2*(u[0]*(W_DRIVE - W0) - u[2]*WC), \
                       1/2*(u[1]*WC + u[3]*(W_DRIVE - W0)), \
                       1/2*(-u[0]*WC - u[2]*(W_DRIVE - W0))])

    return value



def spinor(t, n, u):
# Modified fight hand side of eq A.1 - A.4 in May's nEDM thesis
# using eq 3.28, 3.29 as a basis
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im(b)
    import numpy as np
    x = W_DRIVE*t + PHI
    value = np.array ( [ 1/2*(W0*u[1] + WC*np.cos(x)*u[3]) - WC/2*u[2]*np.sin(x), \
                       1/2*(-W0*u[0] - WC*np.cos(x)*u[2]) - WC/2*u[3]*np.sin(x), \
                       1/2*(-W0*u[3] + WC*np.cos(x)*u[1]) + WC/2*u[0]*np.sin(x), \
                       1/2*(W0*u[2] - WC*np.cos(x)*u[0]) + WC/2*u[1]*np.sin(x) ])

    return value

if ( __name__ == '__main__' ):
    main()
