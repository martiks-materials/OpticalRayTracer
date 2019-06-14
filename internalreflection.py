"""Script - internalreflection:
Instantiates a plane SphericalRefraction object from the raytracermodel module
at z = 125mm with zero curvature, but dispersive.

Instantiates a line of rays with initial direction [0.8734988, 0, 1] with
wavelengths varying between 0.42 and 0.70 micrometres. Propagates rays to plane
and depending on final direction, propagates to one of two OutputPlane objects 
at z=0mm and z=500mm.

Plots the traced path of each ray with a colour that corresponds qualitatively 
to its wavelength, done with a RGB tuple whose components are a function of
its wavelength.
"""

import matplotlib.pyplot as plt
import raytracermodel as rtm
import numpy as np

plane1 = rtm.SphericalRefraction([0, 0, 125], 0.00, 1., 1., 10000,
                                 dispersive=True)
out1 = rtm.OutputPlane([0, 0, 500], 1000000) #This receives refracted rays.
out2 = rtm.OutputPlane([0, 0, 0], 1000000) #This receives reflected rays.
for n in range(281):
    rtm.which = 'n1'
    ray = rtm.Ray([n*0.2, 0, 0], [0.8734988, 0, 1], w=0.7 - 0.001*n)
    print "Wavelength: {} micrometres".format(ray.wavelength())
    plane1.propagate_ray(ray)
    if np.dot(ray.k(), np.array([0, 0, 1])) < 0:
        print "Internally reflected"
        out2.propagate_ray(ray)
    else:
        print "Dispersed"
        out1.propagate_ray(ray)
    x = []
    z = []
    for vector in ray.vertices():
        x.append(vector[0])
        z.append(vector[2])
    W = ray.wavelength()
    plt.plot(z, x, color=((0.4-2*W)**2, 0.8-(3.8*W-1.9)**2, 1-(0.4-1.95*W)**2))
    #The RGB colour tuple contains functions of W that were empirically found.
plt.axis([0, 220, -100, 700])
plt.grid(b=None, which='both', color='k', linestyle='-', linewidth='0.4')
plt.plot([125 for n in range(2)], [-100, 700], 'k--')
xy = [50, 350]
s = "Dispersive Glass"
plt.annotate(s, xy, xytext=None, xycoords='data',
             textcoords='data', arrowprops=None)
xy = [145, 0]
s = "Vacuum: \nRefractive index = 1.000"
plt.annotate(s, xy, xytext=None, xycoords='data',
             textcoords='data', arrowprops=None)
plt.xlabel("Distance along optical axis (mm)")
plt.ylabel("Distance in x direction from optical axis (mm)")
plt.title("Total internal reflection in a dispersive medium")
plt.show()
        
    

   