#!/usr/bin/env python

def main():
    import numpy as np
    from circularFringe import larmor
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    ket = np.array([0.7071, 0, 0.7071, 0])
    xProb = []
    yProb = []
    zProb = []
    dt = .001

    for i in range(1,50):
        zProb.append(ket[0]*ket[0] + ket[1]*ket[1])
        xProb.append(1/2 + ket[0]*ket[2] + ket[1]*ket[3])
        yProb.append(1/2 + ket[1]*ket[2] - ket[3]*ket[0])
        ket = larmor(ket, dt, 100, 4)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111, projection='3d')
    plt.title('Odds of measuring spin up along x, y, z')
    ax.scatter(xProb, yProb, zProb)
    ax.set_xlim3d(0,1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(0,1)
    ax.set_xlabel('P(x)')
    ax.set_ylabel('P(y)')
    ax.set_zlabel('P(z)')
    ax.view_init(30,220)

    plt.show()

    return

if ( __name__ == '__main__' ):
    main()
