#!/usr/bin/env python
#
# Draws ramsey fringes for a neutron in an nEDM experiment,
# where the pi/2 pulse is a linearly
# oscillating magnetic field in the x-y plane
#
# Douglas Wong 4/7/17

# Time parameters
PULSE_1_TIME = 4.286        # [seconds]
PULSE_2_TIME = 4.286        # [seconds]
RK_STEP = 0.001       # [seconds]
PRECESS_TIME = 180      # [seconds]

# Some initial parameters
W_STEP = 0.0005    #[rad s^-1]    Step value of w to make ramsey fringes
W_VAL = 183.2   #[rad s^-1]    What w to start with
W_MAX  = 183.3    #[rad s^-1]    What w to end with

W0_VAL = 183.247172  #[rad s^-1]    Static field strength
WL_VAL = 0.732988688   #[rad s^-1]   Oscillating field strength
PHI_VAL_1 = 0  #[rad]          RF pulse inital phase for first pulse

def main():
    from tqdm import tqdm           # For progress bar, see read me
    from decimal import Decimal     # For error handling
    import numpy as np
    import matplotlib.pyplot as plt

    t0 = 0
    n = 4

    wRange = np.arange(W_VAL, W_MAX, W_STEP)
    zProb = []

    ket = np.zeros ( n )

    if ( (Decimal(str(PULSE_1_TIME)) % Decimal(str(RK_STEP)) != 0) \
        or (Decimal(str(PULSE_2_TIME)) % Decimal(str(RK_STEP)) != 0)):
        print("Error: Pulse time resolution too small for RK integrator")
        return

    for wTemp in tqdm(wRange):
        ket[0] = 1       #neutron starts spin up (ket[0] = Re[a0])
        ket = spinPulse(ket, RK_STEP, PULSE_1_TIME, n, wTemp, W0_VAL, WL_VAL, PHI_VAL_1)
        ket = larmor(ket, PRECESS_TIME, W0_VAL, n)

        #spinPulse 2 has to stay in phase with spinPulse 1 while the larmor precession occurs
        phiVal2 = wTemp*PULSE_1_TIME +PHI_VAL_1 + wTemp*PRECESS_TIME
        ket = spinPulse(ket, RK_STEP, PULSE_2_TIME, n,wTemp, W0_VAL, WL_VAL, phiVal2)
        zProb.append(ket[0]*ket[0] + ket[1]*ket[1])
        ket[:] = 0    # Reset ket for next loop

    #Bloch Siegert Shift
    print("The resonant freq is\t", wRange[zProb.index( min(zProb) )], "rad/s" )

    # Plot Stuff
    pentrackW = np.arange(183.2, 183.3, 0.001)
    pentrackZ = np.array([7.818399655348102772e-01,6.534901541871088737e-01,5.028516202876234242e-01,3.351043227707187766e-01,1.559738404445508986e-01,-2.845688640481612228e-02,-2.119379678260573385e-01,-3.882441566351519935e-01,-5.513679796122282761e-01,-6.957140134273124277e-01,-8.162912165345551641e-01,-9.088966035928353060e-01,-9.702801308791237922e-01,-9.982788204225419015e-01,-9.919085919798140694e-01,-9.514052718641996531e-01,-8.782111529020766305e-01,-7.749092043577318067e-01,-6.451121786223584786e-01,-4.933172330864776800e-01,-3.247376216608716115e-01,-1.451215309315782576e-01,3.943499485539350702e-02,2.226776843790268523e-01,3.983856680354193602e-01,5.605654921946128244e-01,7.036467897090242785e-01,8.226758955349825486e-01,9.134999983524854095e-01,9.729311411765270146e-01,9.988778751934360711e-01,9.904332260454702386e-01,9.479110033493275411e-01,8.728276165385905339e-01,7.678322918654136853e-01,6.365935114396101824e-01,4.836525074671476632e-01,3.142551543818954407e-01,1.341718087046287200e-01,-5.048865325422194350e-02,-2.334678608887226503e-01,-4.085498084793193208e-01,-5.697559328273176238e-01,-7.115422298949749536e-01,-8.289940926992962478e-01,-9.180107291104983025e-01,-9.754680102504504902e-01,-9.993473688897007712e-01,-9.888196910372620829e-01,-9.442768436268774712e-01,-8.673088138652780277e-01,-7.606301234083969742e-01,-6.279638567110044090e-01,-4.738942714001040390e-01,-3.036990639555826599e-01,-1.231702444003499042e-01,6.157079645101940640e-02,2.442618831511327859e-01,4.186923605457313902e-01,5.788993718824667623e-01,7.193664789271760895e-01,8.352196188286595824e-01,9.224115917172838186e-01,9.778834747376506931e-01,9.996905070578797625e-01,9.870817304305592454e-01,9.405266947547776812e-01,8.616880251162483129e-01,7.533442501499623134e-01,6.192716918166167872e-01,4.640964299022360962e-01,2.931270844504537121e-01,1.121767500095643688e-01,-7.262102450170999923e-02,-2.550005877131000243e-01,-4.287571344456804212e-01,-5.879443156568555651e-01,-7.270743894888570535e-01,-8.413151720180345183e-01,-9.266743764621633295e-01,-9.801593201853561199e-01,-9.998995564996788765e-01,-9.852221206489639727e-01,-9.366734141814706538e-01,-8.559873428029171061e-01,-7.460048069891913602e-01,-6.105537312583686482e-01,-4.543006262108660875e-01,-2.825840128876320745e-01,-1.012374195847456498e-01,8.359385206778342792e-02,2.656410304878425999e-01,4.387056553680682702e-01,5.968586527388440865e-01,7.346420009744317747e-01,8.472665614895279473e-01,9.307960345788969247e-01,9.823047018202689662e-01,9.999965840337161627e-01,9.832761838592067871e-01,9.327655812837349902e-01])

    pentrackZ =  (pentrackZ + 1)/2

    plt.plot(wRange, zProb, color = 'black')
    plt.plot(pentrackW, pentrackZ, '.', color = 'cyan')

    plt.grid(True)
    plt.xlabel('w [Hz]')
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
