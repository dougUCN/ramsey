Some random programs hopefully to better understand the ramsey cycle in nEDM experiments
Requires the tqdm library for a progress bar (https://github.com/tqdm/tqdm)

rkCircular.py -- Applies a circularly oscillating B field on the x-y plane
                      and a uniform B field along z to a neutron ket state
rkLinear.py -- Applies a linearly oscillating B field in x-y and a uniform B
                  field along z to a neutron ket state
circularFringe.py -- Generates ramsey fringes using a circular RF field
linearFringe.py -- Generates ramsey fringes using a linear RF field
larmorTest.py -- Visualizes larmor precession given a ket state
blochSiegert.py -- Takes a linear optimized ramsey fringes and finds the Bloch
                    Siegert shift as a function of pulse time
phaseAngle.py -- Plots Bloch Siegert shift as a function of initial phase
angle
asymmetricPulse.py -- Illustrates shift in res. freq when pi/2 pulses are asymmetrical
rkCircularAnimate.py -- Creates an animation of ramsey sequence (circular RF)

circularFringe and linearFringe run slowly, and blochSiegert takes eons.
This is because I'm an idiot and should have used C++ instead,
and to compensate I use a progress bar
(https://github.com/tqdm/tqdm)

Runge Kutta integrator (rk4.py) taken from
https://people.sc.fsu.edu/~jburkardt/cpp_src/rk4/rk4.html
Used only in rkCircular.py and rkLinear.py

RF parameters.txt is a list of values utilized in the LANL nEDM experiment

Many of my comments refer to Daniel May's 1999 Thesis on the ILL nEDM
https://inis.iaea.org/search/searchsinglerecord.aspx?recordsFor=SingleRecord&RN=30045807
Email me for a copy
