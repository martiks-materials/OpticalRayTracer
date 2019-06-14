"""Script - lensoptimize:
Using the curvatures of setup 2 of the planoconvex lenses as an initial guess 
for the optimiser and bounds of -0.1/mm and 0.1/mm for both curvatures, 
minimises the RMS at the paraxial focus of setup 2's planoconvex lens by
varying the curvatures of the lens such that the RMS at the paraxial 
focus is of lowest value.

The image distance is the distance from the axial intersection of the
input-facing surface to the paraxial focus point along the optical (z) axis.

For this image distance, prints the optimum curvatures as a list, with the 
first element being the curvature of the input-facing surface and the second 
element being the curvature of the output-facing surface.

Then prints the minimised RMS deviation from the optical axis at the paraxial
focus of the planoconvex setup 2 when the optimum curvatures have been used.

Optimiser used: scipy.optimize.fmin_tnc:
            -- Function optimised, func=RMS
            -- Initial guess, x0=[0.02, 0.00]
            -- Parameter boundaries, bounds=[(-0.1, 0.1), (-0.1, 0.1)]
            -- Maximum function evaluations, maxfun=100000
            -- Function tolerance, ftol=0.000000001
            -- Approximates function gradient, approx_grad=True
            -- Extra arguments for image distance of 98.45mm, args=(image_dist)
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import raytracermodel as rtm
import numpy as np
import scipy.optimize as spo
import genpolar

def RMS(curv_list, image_dist):
    """For a singlet lens with given surface curvatures, return RMS deviation 
    from optical axis at a given image distance of a uniformly distributed 
    beam of rays of initial diameter 10mm.The beam propagates through the
    singlet lens with refractive index 1.5168 and radius aperture 50mm.
    
    Parameters: curv_list -- List of curvatures for both side of the singlet 
                             lens. The first element of the list is assigned 
                             as the curvature of the surface facing the input,
                             and the second element is assigned as the
                             curvature of the surface facing the output.
               image_dist -- Float image distance from the axial intersection
                             of the input-facing surface to the point of 
                             desired focus.
    """
    curv1 = curv_list[0]
    curv2 = curv_list[1]
    surface1 = rtm.SphericalRefraction([0., 0., 100.], curv1, 1., 1.5168, 50.)
    surface2 = rtm.SphericalRefraction([0., 0., 105.], curv2, 1.5168, 1., 50.)
    output = rtm.OutputPlane([0., 0., image_dist + 100.], 10000.)
    sq_deviation = []
    for r, theta in genpolar.rtuniform(5, 5, 5):
        ray = rtm.Ray([r*np.cos(theta), r*np.sin(theta), 0.], [0., 0., 1.])
        surface1.propagate_ray(ray)
        surface2.propagate_ray(ray)
        output.propagate_ray(ray)
        if ray._terminate == False:
            r_squared = ray.p()[0]**2. + ray.p()[1]**2.
            sq_deviation.append(r_squared)
        else:
            print "RAY TERMINATED\nRadius: {}mm\nAngle: {} ".format(r, theta)
            return None   
    return np.sqrt(np.average(sq_deviation))
    

image_dist = 98.452812504
guess = [0.02, 0.00]
boundaries = [(-0.1, 0.1), (-0.1, 0.1)]
z = spo.fmin_tnc(RMS, x0=guess, args=([image_dist]), approx_grad=True, 
                 bounds=boundaries,  maxfun=100000, ftol=0.000000001)
print"Optimum curvatures:\nSurface 1: {}/mm\nSurface 2: {}/mm".format(z[0][0], 
                                                                       z[0][1])
print "Minimised RMS: {} mm".format(RMS(z[0], image_dist))