#!/usr/bin/env python
#
# Draws ramsey fringes for a neutron in an nEDM experiment,
# where the pi/2 pulse is a circularly
# rotating magnetic field in the x-y plane
#
# Douglas Wong 4/7/17

# Time parameters
PULSE_1_TIME = 3        # [seconds]
PULSE_2_TIME = 3        # [seconds]
RK_STEP = 0.001       # [seconds]
PRECESS_TIME = 10      # [seconds]

# Some initial parameters
W_STEP = 0.01    #[rad s^-1]    Step value of w to make ramsey fringes
W_VAL = 186   #[rad s^-1]    What w to start with
W_MAX  = 190      #[rad s^-1]    What w to end with

W0_VAL = 188  #[rad s^-1]    Static field strength
WC_VAL = 0.5   #[rad s^-1]   Rotating field strength
PHI_VAL_1 = 0  #[rad]        RF pulse inital phase for first pulse

def main():
    from tqdm import tqdm       # For progress bar, see read me
    from decimal import Decimal # For error handling
    import numpy as np
    import matplotlib.pyplot as plt
    t0 = 0
    n = 4

    wRange = np.arange(W_VAL, W_MAX + W_STEP, W_STEP)
    zProb = []

    ket = np.zeros ( n )

    if ( (Decimal(str(PULSE_1_TIME)) % Decimal(str(RK_STEP)) != 0) \
        or (Decimal(str(PULSE_2_TIME)) % Decimal(str(RK_STEP)) != 0)):
        print("Error: Pulse time step too small for RK integrator")
        return

    for wTemp in tqdm(wRange):
        ket[0] = 1       # neutron starts spin up (ket[0] = Re[a0])
        ket = spinPulse(ket, RK_STEP, PULSE_1_TIME, n, wTemp, W0_VAL, WC_VAL, PHI_VAL_1)
        ket = larmor(ket, PRECESS_TIME, W0_VAL, n)

        #spinPulse 2 has to stay in phase with spinPulse 1 while the larmor precession occurs
        phiVal2 = wTemp*(PULSE_1_TIME) + PHI_VAL_1 + wTemp*PRECESS_TIME
        ket = spinPulse(ket, RK_STEP, PULSE_2_TIME, n, wTemp, W0_VAL, WC_VAL, phiVal2)
        zProb.append(ket[0]*ket[0] + ket[1]*ket[1])
        ket[:] = 0    # Reset ket for next loop

    print("The resonant freq is\t", wRange[zProb.index( min(zProb) )], " rad/s" )

    # Plot Stuff
    plt.plot(wRange,zProb)
    plt.xlabel('w [rad s^-1]')
    plt.ylabel('P(z)')
    plt.show()

    return

def spinPulse(u0, dt, tmax, n, w, w0, wc, phi):
# Runge Kutta integration for SPINOR with step size dt
# Requires a [n] element vector u0 with inital parameters
# Returns the updated array u0 specifying where the neutron has ended up
    import numpy as np
    tRange = np.arange(0, tmax, dt)
    for tTemp in tRange:
        u0 = rkStep ( tTemp, n, u0, dt, w, w0, wc, phi, spinor)
    return u0

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


def spinor(t, n, u, w, w0, wc, phi):
# Modified right hand side of eq A.1 - A.4 in May's nEDM thesis
# using eq 3.28, 3.29 as a basis
# u[0] = Re(a'), u[1] = Im(a'), u[2] = Re(b'), u(3) = Im(b')
    import numpy as np
    if ( n != 4 ):
        return np.zeros(n)

    x = w*t + phi
    value = np.array ( [ 1/2*(w0*u[1] + wc*np.cos(x)*u[3]) - wc/2*u[2]*np.sin(x), \
                       1/2*(-w0*u[0] - wc*np.cos(x)*u[2])  - wc/2*u[3]*np.sin(x), \
                       1/2*(-w0*u[3] + wc*np.cos(x)*u[1]) + wc/2*u[0]*np.sin(x), \
                       1/2*(w0*u[2] - wc*np.cos(x)*u[0]) + wc/2*u[1]*np.sin(x) ])

    return value

def rkStep ( t0, m, u0, dt, w, w0, wc, phi, f):
# Takes one runge kutta step
# adapted from rk4 library to accomodate an f with frequency and phase

  import numpy as np
#
#  Get four sample values of the derivative.
#
  f0 = f ( t0, m, u0, w, w0, wc, phi )

  t1 = t0 + dt / 2.0
  u1 = np.zeros ( m )
  u1[0:m] = u0[0:m] + dt * f0[0:m] / 2.0
  f1 = f ( t1, m, u1, w, w0, wc, phi )

  t2 = t0 + dt / 2.0
  u2 = np.zeros ( m )
  u2[0:m] = u0[0:m] + dt * f1[0:m] / 2.0
  f2 = f ( t2, m, u2, w, w0, wc, phi )

  t3 = t0 + dt
  u3 = np.zeros ( m )
  u3[0:m] = u0[0:m] + dt * f2[0:m]
  f3 = f ( t3, m, u3, w, w0, wc, phi )
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
