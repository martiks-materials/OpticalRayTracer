"""Script - beamtracing:
Instantiates a spherical refracting surface of curvature +0.02/mm placed a
distance z = 100mm along the optical axis, with refractive indices
1.00 for z < 100mm and 1.5 for z > 100mm. Estimates the paraxial focus 
ofthe surface by propagating a beam 10e-4mm from optical axis and parallel
to it through the surface and to an OutputPlane object at z = 250mm. Moves
the OutputPlane to this point and propagates a uniform beam from 
the input plane at z = 0mm through the surface to this output.

Plots spot diagrams at the input plane and the output plane.

Additionally, plots the path of each ray in the beam in a 3D diagram.
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import raytracermodel as rtm
import genpolar
import numpy as np


sphere = rtm.SphericalRefraction([0, 0, 100], 0.03, 1, 1.5, 33)
test_output = rtm.OutputPlane([0, 0, 250], 33)
test_ray = rtm.Ray([0.0001, 0, 0], [0, 0, 1])
sphere.propagate_ray(test_ray)
test_output.propagate_ray(test_ray)
p1 = test_ray.vertices()[1]
p2 = test_ray.vertices()[2]
gradient =((p1[0]-p2[0])/(p1[2]-p2[2]))
#This is the gradient of the line equation for the ray's final path.
paraxial_focus = p1[2] - p1[0]/gradient
print "Paraxial Focus:  {} mm ".format(paraxial_focus)
output = rtm.OutputPlane([0, 0, paraxial_focus], 33)
x1 = []
y1 = []
x2 = []
y2 = []
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
for r, t in genpolar.rtuniform(10, 10, 10):                      
    x0 = []
    y0 = []
    z0 = []
    ray = rtm.Ray([r*np.cos(t), r*np.sin(t), 0], [0, 0, 1])
    x1.append(ray.p()[0])
    y1.append(ray.p()[1])
    sphere.propagate_ray(ray)
    output.propagate_ray(ray)
    for vector in ray.vertices():
        x0.append(vector[0])
        y0.append(vector[1])
        z0.append(vector[2])
    ax.plot(z0, x0, y0, 'ro--')
    x2.append(ray.p()[0])
    y2.append(ray.p()[1])
plt.title("3D plot of Traced rays")
ax.set_xlabel("Distance along optical axis (mm)")
ax.set_ylabel("x distance from optical axis (mm)")
ax.set_zlabel("y distance from optical axis (mm)")
plt.figure(2)
plt.plot(x1, y1, 'go')
plt.grid(b=None, which='both', color='r', linestyle='-', linewidth='0.4') 
plt.title("Spot Diagram at Input Plane z = 0")
plt.xlabel("X")
plt.ylabel("Y")
plt.figure(3)
plt.plot(x2, y2, 'go')
plt.grid(b=None, which='both', color='r', linestyle='-', linewidth='0.4')
plt.title("Spot Diagram at Output Plane z = 200")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

sq_deviation = []
for i in range(len(x2)): 
    sq_deviation.append(x2[i]**2 + y2[i]**2)
rms_deviation = np.sqrt(np.average(sq_deviation))
#This determines the RMS deviation of all the rays in the uniform beam.
print "Root Mean Square Deviation: {} mm ".format(rms_deviation)

    

    


