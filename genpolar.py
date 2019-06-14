"""genpolar module:
Generate positions in polar coordinates uniformly distributed
within a circle of specified radius. Useful for producing initial
positions for a collimated beam of Ray objects in the ray
tracer model.
"""

import numpy as np

def rtpairs(R, N):
    """Generate a sequence of pairs of radii and angles.
    
    Parameters: R -- List of radii r.
                N -- List of number of angles for each r.
                     N must have the same length as R.
    
    Raises: Exception -- If lists aren't of the same length.
                
    Yield pairs of (radius, angle) upon iteration.
    """
    if len(N) != len(R):
        raise Exception("Input lists must have the same length.")
    for i in range(len(R)):
        no_of_angles = N[i]
        spacing = 2*np.pi/no_of_angles
        angles = np.arange(0, 2*np.pi, spacing)
        for t in angles:
            r = R[i]
            yield r, t
            
def rtuniform(n, rmax, m):
    """Generate a sequence of (radius, angle) pairs that
    are uniformly distributed over a circle of given radius.
    
    Parameters: n -- Int number of rings required.
             rmax -- Maxiumum radial extent of points.
                m -- Int number such that the ith ring
                     from the centre contains m*i points.
                     
    Yield a sequence of (radius, angle) pairs.
    """
    R = [i*rmax/n for i in range(n+1)]
    N = [1]+[ m*i for i in range(1, n+1)]
    for r,t in rtpairs(R, N):
        yield r, t
        


