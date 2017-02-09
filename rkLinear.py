#! /usr/bin/env python
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
    tmax = 2*np.pi
    t0 = 0

    u0 = np.zeros ( n )
    u0[0] = 1
    u0[1] = 0
    u0[2] = 0
    u0[3] = 0

    time = [0]
    aNorm = []

    i = 0
    while ( True ):
        # print ( '  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g' \
        #         % ( i, time[i], u0[0], u0[1], u0[2], u0[3]))

        # odds of measuring spin up along z
        aNorm.append( np.sqrt(np.power(u0[0], 2) + np.power(u0[1], 2)) )
        # aNorm.append(u0[1])
  #
  #  Stop if we've exceeded TMAX.
  #
        if ( tmax <= t0 ):
            break

        i = i + 1
  #
  #  Otherwise, advance to time T1, and have RK4 estimate
  #  the solution U1 there.
  #
        t1 = t0 + dt
        u1 = rk4vec ( t0, n, u0, dt, spinor)     #takes one RK step
  #
  #  Shift the data to prepare for another step.
  #
        t0 = t1
        u0 = u1.copy ( )
        time.append(t1)
  #  End of While loop

  #  Graphing stuff
    plt.plot(time, aNorm)
    plt.title('Odds of measuring spin up along z')
    plt.xlabel('time [s]')
    plt.ylabel('P(z)')
    plt.axis([0,2*np.pi,-1,1.5])
    plt.show()
    print ( '' )
    print ( 'RKLINEAR:' )
    print ( '  Normal end of execution.' )
    return

def spinor(t, n, u):
# Right hand side of eq A.1 - A.4 in May nEDM thesis
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im b
    import numpy as np
    w = 188   #[rad s^-1]
    w0 = 188  #[rad s^-1]
    wl = 1.57 #[rad s^-1]
    # wl = 0    #[rad s^-1]
    phi = 0   #[rad]
    value = np.array ( [ 1/2*(w0*u[1] + wl*np.cos(w*t + phi)*u[3]), \
                       1/2*(-w0*u[0] - wl*np.cos(w*t + phi)*u[2]), \
                       1/2*(-w0*u[3] + wl*np.cos(w*t + phi)*u[1]), \
                       1/2*(w0*u[2] + wl*np.cos(w*t + phi)*u[0])])

    return value

def testFunc(t, n, u):
# Driven, damped pendulum. For testing. u[0]= theta, u[1] = w
    import numpy as np
    # value = np.array ( [ u[1], \
    #                    -0.5*u[1]-np.sin(u[0])+0.5*np.cos(0.7*t)])
    value = np.array ([1/2*188*u[1], 1/2*(-188)*u[0]])
    return value

if ( __name__ == '__main__' ):
    main()
