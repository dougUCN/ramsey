#!/usr/bin/env python
#
# What happens when the pi/2 pulses are assymetric
#
# Douglas Wong 5/28/17

# Time parameters
PULSE_2_INIT = 4.286        # [seconds]
PULSE_2_FINAL = 4.286     #[seconds]
PULSE_STEP = 0.001    # [rad] Time step for x axis of final graph
RK_STEP = 0.001       # [seconds] For Runge Kutta integrator
PRECESS_TIME = 180      # [seconds]

# Some initial parameters
W_STEP = 0.0000005    #[rad s^-1]    Step length of search around w0
W_STEP_NUM = 100    # Number of steps to search around w0
W0_VAL = 183.247  #[rad s^-1]    Static field strength
WL_VAL = None  #[rad s^-1]    RF field strength
PULSE_1 = 4.286  #[seconds]  Time in which pi/2 pulse applied
PHASE_INIT = 0 #[rad]

def main():
    from tqdm import tqdm   # Magic progress bar
    from decimal import Decimal
    import numpy as np
    import matplotlib.pyplot as plt
    t0 = 0
    n = 4

    wRange = np.arange(W0_VAL - W_STEP*W_STEP_NUM, W0_VAL + W_STEP*W_STEP_NUM, W_STEP)
    pulseRange = np.arange(PULSE_2_INIT, PULSE_2_FINAL + PULSE_STEP, PULSE_STEP)
    zProb = []
    rkError = []
    shift = []

    ket = np.zeros ( n )

    # Optimized ramsey resonance calculation
    if WL_VAL is None:
        wl = np.pi / PULSE_1  # Calculate w_l for optimized ramsey resonance
        pulse1 = PULSE_1
    elif PULSE_1 is None:
        pulse1 = np.pi / WL_VAL
        wl = WL_VAL
    else:
        print("Warning: not calculating optimized ramsey fringes")
        pulse1 = PULSE_1
        wl = WL_VAL

    # Error checking
    if (Decimal(str(pulse1)) % Decimal(str(RK_STEP)) != 0):
        print("Error: RK step too large for pulse time")
        print("Pulse time: ", PULSE_1)
        return

    for pulse2 in tqdm(pulseRange):
        for wTemp in wRange:    # Does one optimized ramsey fringe for given pulse time
            ket[0] = 1       #neutron starts spin up (ket[0] = Re[a0])
            ket = spinPulse(ket, RK_STEP, PULSE_1, n, wTemp, W0_VAL, wl, PHASE_INIT)
            ket = larmor(ket, PRECESS_TIME, W0_VAL, n)

            #spinPulse 2 has to stay in phase with spinPulse 1 while the larmor precession occurs
            phiVal2 = wTemp*PULSE_1 + PHASE_INIT + wTemp*PRECESS_TIME
            ket = spinPulse(ket, RK_STEP, pulse2, n,wTemp, W0_VAL, wl, phiVal2)
            zProb.append(ket[0]*ket[0] + ket[1]*ket[1])
            ket[:] = 0    # Reset ket for next loop

        #Resonant frequency
        shift.append(W0_VAL - wRange[zProb.index( min(zProb) )])
        zProb.clear()   #clear for next loop


    # Plot Stuff
    np.savetxt('asymPulse.txt', np.c_[pulseRange - PULSE_1, np.divide(shift, 2*np.pi)])
    plt.plot(pulseRange - PULSE_1, np.divide(shift, 2*np.pi))
    plt.grid(True)
    plt.xlabel('Difference in pulse time [s]')
    plt.ylabel('Shift in res freq [Hz]')
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
# u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im(b)
    import numpy as np
    if ( n != 4 ):
        return np.zeros(n)

    x = w*t + phi
    value = np.array ( [ 1/2*(w0*u[1] + wc*np.cos(x)*u[3]), \
                       1/2*(-w0*u[0] - wc*np.cos(x)*u[2]), \
                       1/2*(-w0*u[3] + wc*np.cos(x)*u[1]), \
                       1/2*(w0*u[2] - wc*np.cos(x)*u[0]) ])

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
