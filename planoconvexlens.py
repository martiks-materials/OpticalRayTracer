"""Script - planoconvexlens:
Defines two setups of OpticalElement derived objects that characterise the two 
orientations of a planoconvex lens with refractive index 1.5168 and one curved 
side of curvature +/-0.02/mm. The code runs both setups without the user 
changing any parameters midway. For each setup, estimates the paraxial focus 
using a parallel ray 10e-4mm from the optical axis and propagates it through
the lens to an OutputPlane object initially at 200mm. The output plane is then 
moved to the paraxial focus and a uniform circular beam of diameter 10mm is 
propagated to the lens.

This script prints the paraxial focus of the lens with respect to the z = 0mm
input plane, and the root mean square deviation from the optical axis at the 
paraxial focus for all the Ray objects in the beam.

Plots spot diagrams for at the paraxial focus plane for both orientations
of the lens.

Reduces beam diameter and runs the propagation again, recording only RMS and
diffraction scale, which are both plotted on a graph against beam diameter 
for both setups.

Finds and prints the maximum beam diameter where the RMS and diffraction scale
are within 5e-5mm of each other. Also prints the RMS that this corresponds to.
"""

import matplotlib.pyplot as plt
import numpy as np
import raytracermodel as rtm
import genpolar

plane1 = rtm.SphericalRefraction([0, 0, 100], 0.0, 1., 1.5168, 150)
convex1 = rtm.SphericalRefraction([0, 0, 105], -0.02, 1.5168, 1., 150)
plane2 = rtm.SphericalRefraction([0, 0, 105], 0.0, 1.5168, 1., 150)
convex2 = rtm.SphericalRefraction([0, 0, 100], 0.02, 1., 1.5168, 150)
output = rtm.OutputPlane([0, 0, 200], 100)

Setup1 = [plane1, convex1, output] #Here, the plane surface faces the input.
Setup2 = [convex2, plane2, output] #Here, the curved surface faces the input.
Experiments = [Setup1, Setup2]

for setup in Experiments:
    print "Setup {}:".format(Experiments.index(setup) + 1)
    ray = rtm.Ray([0.0001, 0, 0], [0, 0, 1]) 
    for element in setup:
        element.propagate_ray(ray)
    p1 = ray.vertices()[2]
    p2 = ray.vertices()[3]
    gradient =((p1[0]-p2[0])/(p1[2]-p2[2])) 
    #Finds line equation of ray's path between lens and the output plane.
    if gradient >= 0:
        print "LENS CANNOT FOCUS"
    paraxial_focus = p1[2] - p1[0]/gradient
    print "Paraxial Focus:  {} mm ".format(paraxial_focus)
    setup[2] = rtm.OutputPlane([0, 0, paraxial_focus], 100)
    #Moves output to paraxial focus.
    sq_deviation = []
    plt.figure(1)
    plt.subplot(1, 2, Experiments.index(setup) + 1)
    for r, t in genpolar.rtuniform(5, 5, 10):
        ray = rtm.Ray([r*np.cos(t), r*np.sin(t), 0], [0, 0, 1])                             
        for element in setup:
            element.propagate_ray(ray)
        if ray._terminate == False:
            plt.plot(ray.p()[0], ray.p()[1], 'ro')
            sq_deviation.append(ray.p()[0]**2+ray.p()[1]**2)
    plt.grid(b=None, which='both', color='c', linestyle='-', linewidth='0.4')
    plt.title("Spot diagram at paraxial focus")
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    plt.axis([-0.08, 0.08, -0.07, 0.09])
    s = "Setup {}\nParaxial Focus {} mm\nRMS: {} mm \nBeam Diameter:\
    10 mm".format(Experiments.index(setup) + 1, "%.5G"%(paraxial_focus),\
                  "%.5G"%np.sqrt(np.average(sq_deviation)))
    xy = [-0.04, 0.06]
    plt.annotate(s, xy, xytext=None, xycoords='data', textcoords='data',
                 arrowprops=None)
    print "RMS Deviation:  {} mm ".format(np.sqrt(np.average(sq_deviation)))
    diameter = []
    RMS = []
    x_limit = []
    for n in np.logspace(-3, 1, 50):
        sq_dev = []
        for r, t in genpolar.rtuniform(10, 0.5*n, 5):
            ray = rtm.Ray([r*np.cos(t), r*np.sin(t), 0], [0, 0, 1])                              
            for element in setup:
                element.propagate_ray(ray)
            if ray._terminate == False:
                sq_dev.append(ray.p()[0]**2+ray.p()[1]**2)
        diameter.append(n)
        RMS.append(np.sqrt(np.average(sq_dev)))
        f = paraxial_focus - 100
        #Focal length f is the focus position minus the lens position.
        x = 600e-9*f/n
        #This approximates the order of magnitude of the diffraction limit.
        x_limit.append(x)
    limited_lens = []
    for index in range(len(RMS)):
        if abs(RMS[index] - x_limit[index]) < 0.00005:
            limited_lens.append(diameter[index])
    limited_diameter = np.max(limited_lens)
    limited_rms = RMS[diameter.index(limited_diameter)]
    print "Approximate beam diameter when lens is diffraction limited: {} mm"\
           .format(limited_diameter)
    print "Approximate diffraction limited RMS: {} mm".format(limited_rms)
    plt.figure(2)
    plt.subplot(1, 2, Experiments.index(setup) + 1)
    plt.plot(diameter, RMS, 'b-')
    plt.plot(diameter, x_limit, 'g-')
    plt.plot([limited_diameter for n in range(2)], [0, 0.0015], 'k--')
    plt.plot([0, 3], [limited_rms for n in range(2)], 'k--')
    plt.grid(b=None, which='both', color='c', linestyle='-', linewidth='0.4')
    plt.axis([0, 3, 0, 0.0015])
    plt.title("Setup {}".format(Experiments.index(setup) + 1))
    plt.xlabel("Beam Diameter (mm)")
    plt.ylabel("RMS/diffraction Scale (mm)")
    plt.legend(['RMS Deviation', 'Diffraction Scale'], loc='upper center')
plt.show()
     
           

