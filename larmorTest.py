#!/usr/bin/env python

def main():
    import numpy as np
    from circularFringe import larmor
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    ket = np.array([-0.82573592,  0.32626352,  0.04390854, -0.4580222 ])
    xProb = []
    yProb = []
    zProb = []
    dt = .01
    w0 = 20

    for i in range(1,200):
        zProb.append(ket[0]*ket[0] + ket[1]*ket[1])
        xProb.append(1/2 + ket[0]*ket[2] + ket[1]*ket[3])
        yProb.append(1/2 + ket[1]*ket[2] - ket[3]*ket[0])
        ket = larmor(ket, dt, w0, 4)

    print(ket)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111, projection='3d')
    ax.plot(xProb, yProb, zProb, color="#005F73")
    ax.set_xlim3d(0,1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(0,1)
    ax.set_xlabel('P(x)')
    ax.set_ylabel('P(y)')
    ax.set_zlabel('P(z)')
    ax.view_init(12,220)

    plt.show()

    return

if ( __name__ == '__main__' ):
    main()
