#!/usr/bin/env python
#
#   4th order runge kutta on linear on nEDM in
#   a Magnetic Field with Uniform Static and Transverse
#   Linearly Oscillating Components
#
#   Douglas Wong 2/8/17


def main():
    from rk4 import rk4vec
    from rk4 import rk4vec_test_f
    import numpy as np
    import matplotlib.pyplot as plt

    # Initialize stuff
    n = 4     # number of equations in the vector
    dt = 0.001
    tmax = 2
    t0 = 0

    u0 = np.zeros ( n )
    u0[0] = 1
    u0[1] = 0
    u0[2] = 0
    u0[3] = 0

    time = [0]
    zProb = []
    xProb = []
    yProb = []

    i = 0
    while ( True ):
        # print ( '  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g' \
        #         % ( i, time[i], u0[0], u0[1], u0[2], u0[3]))

        # odds of measuring spin along z, x, and y
        # u0[0] = Re[a], u0[1] = Im[a], u0[2] = Re[b], u0[3] = Im[b]

        zProb.append(np.power(u0[0], 2) + np.power(u0[1], 2))
        xProb.append(1/2 + u0[0]*u0[2] + u0[1]*u0[3])
        yProb.append(1/2 + u0[1]*u0[2] - u0[3]*u0[0])

        if ( tmax <= t0 ):
            break

        i = i + 1

        # Takes one RK step
        t1 = t0 + dt
        u1 = rk4vec ( t0, n, u0, dt, spinor)     #takes one RK step

        # Shift data
        t0 = t1
        u0 = u1.copy ( )
        time.append(t1)
  #  END WHILE

    plotStuff(xProb, yProb, zProb, time)

    print ( '' )
    print ( 'RKLINEAR:' )
    print ( '  Normal end of execution.' )
    return

def spinor(t, n, u):
# Right hand side of eq A.1 - A.4 in May nEDM thesis
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im b
    import numpy as np
    w = 20   #[rad s^-1]
    w0 = 20  #[rad s^-1]
    wl = 3.1415 #[rad s^-1]
    phi = 0   #[rad]
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
