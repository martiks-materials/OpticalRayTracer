"""Ray Tracer Module:

Provides Ray class which models a travelling ray, along with an optical element
base class from which the SphericalRefraction, SphericalReflection, and 
OutputPlane classes inherit from. This is to model optics experiments with 
refractive surfaces.            
"""

import numpy as np
import numpy.linalg as la

which = 'n1'
    
def unit(vector):
    """Input 3-dimensional numpy array. Return the normalised form of the
    vector as a 3-dimensional numpy array.
    """
    if la.norm(vector) == 0:  #For the zero vector, return the zero vector. 
        return vector
    else: return vector*la.norm(vector)**(-1.)
    
def reflect(incident, normal):
    """Reflect a given incident ray from a given surface normal.
    
    Parameters: incident -- 3-dimensional numpy array of incident direction.
                  normal -- 3-dimensional numpy array of surface normal.
             
    Return 3-dimensional numpy array representing the new direction vector.
    """
    norm = unit(np.array(normal))
    v_inc = unit(np.array(incident))
    costheta = np.dot(norm, v_inc)
    v_rfl = v_inc - (2*costheta)*norm
    return v_rfl

def snells(incident, normal, n1, n2, dtype=float):
    """Return the refracted direction of a ray upon contact with an interface 
    between regions of differing refractive indices.
    
    Parameters: incident -- 3-dimensional numpy array of incident direction.
                  normal -- 3-dimensional numpy array of surface normal.
                      n1 -- Float refractive index where incident ray leaves.
                      n2 -- Float refractive index where ray is refracting.
    
    Return new direction of ray in the form of a 3-dimensional numpy array.
    """
    norm = unit(np.array(normal))
    v_inc = unit(np.array(incident))
    costheta_1 = np.dot(norm, -v_inc)
    sintheta_1 = np.sqrt(abs(1.0 - (costheta_1)**2))
    if sintheta_1 > (n2/n1): 
        #If the angle is beyond critical, the ray will internally reflect.
        return reflect(v_inc, norm)
    theta_2 = np.arcsin(n1/n2*sintheta_1)
    if sintheta_1 == 0.:    #For zero incident angle, direction is unchanged.
        return v_inc 
    if costheta_1 > 0.:
        v_out = (n1/n2)*v_inc+ ((n1/n2)*costheta_1 - np.cos(theta_2))*norm     
    else:              #Form of v_out depends on direction of surface normal.
        v_out = (n1/n2)*v_inc+ ((n1/n2)*costheta_1 + np.cos(theta_2))*norm
    return v_out

def index(W, A=[1.03961212, 0.231792344, 1.01046945], 
          B=[6.00069867e-3, 2.00179144e-2, 1.03560653e2]):
    """For a particular wavelength, return materials effective refractive index
    as a float.
    
    Parameters: W -- Float wavelength quoted in micrometres.
                A -- List of material Sellmeier coefficients A = [A1, A2, A3]
                     (Default A=[1.03961212, 0.231792344, 1.01046945]).
                B -- List of material Sellmeier coefficients B = [B1, B2, B3]
                     (Default B=[6.00069867e-3, 2.00179144e-2, 1.03560653e2]).
                     
    
    Default coefficients represent those of borosilicate crown glass.
    """
    n = 1.
    for i in range(3):
        n += (A[i]*W**2)/(W**2 - B[i])
    index = np.sqrt(n)
    return index
            
    
class Ray:
    """Class for modelling a travelling ray:
    
    Initialization: p -- Initial point as 3-dimensional numpy array.
                    k -- Initial direction as 3-dimensional numpy array
           **kwargs w -- wavelength in micrometres float (Default 0.58).
                                         
    Raises: Exception -- If p input does not have 3 components.
                      -- If k input does not have 3 components.                                   
                                        
    Other properties: -- List of points where the ray has changed directions.
                      -- Termination (Boolean): If the ray cannot be propagated
                         then it is considered terminated (Default False).
    """
    def __init__(self, p=[0., 0., 0.], k=[0., 0., 0.], w=0.58):
        if len(p) != 3:
            raise Exception("p must have 3 components.")
        if len(k) != 3:
            raise Exception("k must have 3 components.")
        self.__p = np.array(p, dtype=float)
        self.__k = np.array(k, dtype=float)
        self.__w = w
        self.__points = []
        self.__points.append(self.__p)
        self._terminate = False
    
    def p(self):
        """Return the current point of the ray in space as a
        3-dimensional numpy array.
        """
        return self.__p
        
    def k(self):
        """Return the current direction of the ray as a 
        3-dimensional numpy array.
        """
        return self.__k
    
    def wavelength(self):
        """Return the wavelength property of the ray as a float."""
        return self.__w
        
    def append(self, p, k=None):
        """Update the current position and direction, adding the new position
        to the internal list of points of the ray.
        
        Parameters: p -- New position as 3-dimensional numpy array.
                    k -- New direction as 3-dimensional numpy array.
        
        """
        self.__p = np.array(p, dtype=float)
        if k != None:
            self.__k = np.array(k, dtype=float)
        self.__points.append(self.__p)
        
    def __repr__(self):
        return "Ray({}, {})".format(self.__p, self.__k)
        
    def __str__(self):
        return "Ray:\nDirection: {}\nPosition : {}".format(self.__k, self.__p)

    def vertices(self):
        """Return internal list of visited points of the ray."""
        return self.__points
    
    def terminate(self):
        """Terminate ray, signifying that it has not propagated."""
        self._terminate = True
        
        

class OpticalElement:
    """
    Base class for the experimental components of the model.
    
    Initialization: z -- Position vector where the element intersects
                         the optical axis (Default [0., 0., 0.]).
                   Ra -- Radius of the aperture is the extent of the lens
                         from the optical axis (Default 50.).
                         
    Raises: Exception -- If z input does not have 3 components.
    """
    def __init__(self, z=[0., 0., 0.], Ra=50.):
        if len(z) != 3:
            raise Exception("z must have 3 components.")
        self._z = np.array(z)
        self._Ra = Ra
        
    def propagate_ray(self, ray):
        """Propagate a ray through the element.
           
        Parameters: ray -- Ray object.
           
        Raises: NotImplementedError.
        """
        raise NotImplementedError()                                                         
        

class SphericalRefraction(OpticalElement):
    """Model of an spherical OpticalElement surface which refracts and 
    propagates incoming rays.
    
    Initialization: z -- Position vector as a 3-dimensional numpy array 
                         representing where the element intersects the optical 
                         axis (Default [0., 0., 0.]).
            curvature -- Float representing the extent to which the lens is 
                         convex or concave. Inverse of the radius of curvature.
                         Zero curvature indicates a plane lens(Default 1.).
                   n1 -- Float refractive index on left hand side of surface 
                         (Default 1.).
                   n2 -- Float refractive index on right hand side of surface
                         (Default 1.).
                   Ra -- Float radius of the aperture is the extent of the lens
                         from the optical axis (Default 50.).
                         
    Raises: Exception -- If z input does not have 3 components.
    """
    def __init__(self, z=[0., 0., 0.], curvature=1., n1=1., n2=1., Ra=50.,
                 dispersive=False):
        self.__curvature = float(curvature)
        self.__n1 = n1
        self.__n2 = n2
        self.__dispersive = dispersive
        OpticalElement.__init__(self, z, Ra)
        
    def centre(self):
        """Return the centre of curvature of the surface as a 3-dimensional 
        numpy array. If curvature is zero, return the axial intersection. 
        """
        if self.__curvature == 0.:
                return self._z
                #Defines plane's "centre of curvature" as its axis intesection.
        Rc = (self.__curvature)**(-1.)
        u = np.array([0., 0., 1.])
        return self._z + Rc*u
    
    def disperse_index(self, ray, which='n1'):
        """Change the refractive index either side of a surface based on the
        wavelength of the given ray. Specify the correct
        refractive index to alter by inputting which equal to either 'n1' or
        'n2'.
        
        Parameters: ray -- Ray object.
                  which -- String of either 'n1' or 'n2' (Default 'n1').
        """
        if which is 'n1':
            self.__n1 = index(ray.wavelength())
        elif which is 'n2':
            self.__n2 = index(ray.wavelength())
    
    def intercept(self, ray):
        """Return the point on the surface as a 3-dimensional numpy array 
        where a given ray is intercepted by the SphericalRefraction object.
        
        Parameters: ray -- Ray object.
        """
        length = 0.
        if ray._terminate is True:
            return None
        if self.__curvature == 0.: #This considers the case of a planar lens.
            shift_p = ray.p() - np.array([ray.p()[0], 0., 0.])\
                      - np.array([0., ray.p()[1], 0.])
            theta = np.arccos(np.dot(np.array([0., 0., 1.]), unit(ray.k())))   
            r = shift_p-self._z                                               
            length = abs(la.norm(r)*(np.cos(theta))**(-1.))
            intercept = ray.p() + length*unit(ray.k())
            if np.sqrt(intercept[0]**2. + intercept[1]**2.) > self._Ra:
                #Considers if the ray's intercept exceeds the aperture radius.        
                return None                                                 
            else:
                return intercept
        else: 
            Rc = self.__curvature**(-1.)
            r = ray.p() - self.centre()
        r_dot_k = np.dot(r, unit(ray.k()))
        if (r_dot_k)**2 - np.dot(r, r) + Rc**2 < 0.:
            #If the quadratic "discriminant" is negative, there is no solution.                          
            return None                                                       
        length_1 = -r_dot_k + ((r_dot_k)**2. - np.dot(r, r) + Rc**2.)**(0.5)
        length_2 = -r_dot_k - ((r_dot_k)**2. - np.dot(r, r) + Rc**2.)**(0.5)
        if self.__curvature > 0.:                                               
            length = min(abs(length_1), abs(length_2))
        elif self.__curvature < 0.:                                          
            if la.norm(ray.p() - self.centre()) < abs(Rc):
                #Considers the case where the ray begins inside the sphere.
                if np.dot(ray.p() - self.centre(),\
                          np.array([0., 0., 1.])) >= 0.:
                    length = min(abs(length_1), abs(length_2))
                else:
                    length = max(abs(length_1), abs(length_2))
            else:                                                              
                length = max(abs(length_1), abs(length_2))
        intercept = ray.p() + length*unit(ray.k())
        if np.sqrt(intercept[0]**2. + intercept[1]**2.) > self._Ra:
            #Considers if the ray's intercept exceeds the aperture radius.         
            return None                                                     
        else:                                                                   
            return intercept
            
    def propagate_ray(self, ray):
        """Propagate a given ray that has been intercepted by the surface
        and change the direction of the ray depending on the refractive
        properties on either side of the surface.
        
        Parameters: ray -- Ray object.
        """
        if self.intercept(ray) is None:                                      
            ray.terminate()
            print "Ray did not propagate."
            return None
        if self.__curvature == 0.0:    #Considers case of a planar surface.
            if self._z is np.array([0., 0., 0.]):
                surface_normal = np.array([0, 0, -1])
            else:
                surface_normal = unit(self._z)
        else: 
            surface_normal = unit(self.intercept(ray) - self.centre())
        if self.__dispersive is True:
            #This considers if the SphericalRefraction object is dispersive.                                          
            self.disperse_index(ray, which=which)
            #Updates the appropriate refractive index.                              
        new_direction = snells(ray.k(), surface_normal, self.__n1, self.__n2)
        ray.append(self.intercept(ray), new_direction)
        
        
class OutputPlane(OpticalElement):
    """Model for an OpticalElement that receives and propagates 
    incoming rays incident on it.
    
    Initialization: z -- Intersection of the element with the optical axis as a
                         3-dimensional numpy array (Default [0., 0., 0.]).
                   Ra -- Float radius of the aperture, extent of the element 
                         from the optical axis (Default 50.).
    
    Raises: Exception -- If z input does not have 3 components.
    """
    def intercept(self, ray):
        """Return the point on the surface as a 3-dimensional numpy array where
        a given ray is intercepted by the OutputPlane object.
        
        Parameters: ray -- Ray object.
        """
        if ray._terminate is True:
            return None
        shift_p = ray.p() - np.array([ray.p()[0], 0., 0.])\
                  - np.array([0., ray.p()[1], 0.])
        theta = np.arccos(np.dot(np.array([0., 0., 1.]), unit(ray.k())))
        r = shift_p - self._z
        length = abs(la.norm(r)*(np.cos(theta))**(-1.))
        intercept = ray.p() + length*unit(ray.k())
        if np.sqrt(intercept[0]**2. + intercept[1]**2.) > self._Ra:
            #Considers if the ray's intercept exceeds the aperture radius.
            ray.terminate()
            return None 
        else:
            return intercept
            
    def propagate_ray(self, ray):
        """Propagate a given ray that has been intercepted by the surface, and 
        append the new position to the internal vertex list of the ray.
        Does not change direction of the ray.
        
        Parameters: ray -- Ray object.
        """
        if ray._terminate is True:
            return None
        if self.intercept(ray) is None:
            ray.terminate()
            return None
        ray.append(self.intercept(ray))
        
class SphericalReflection(OpticalElement):
    """Model for spherical reflecting optical element with no refractive or
    dispersive properties.
    
    Initialization: z -- Intersection of the element with the optical axis
                         represented by a 3-dimensional numpy array
                         (Default [0., 0., 0.]).
            curvature -- Float representing the extent to which the mirror is 
                         convex or concave. Inverse of the radius of curvature. 
                         Zero curvature corresponds to a plane mirror
                         (Default 1).
                   Ra -- Float radius of the aperture, extent of the element
                         from the optical axis (Default 50.).
                         
    Raises: Exception -- If z input does not have 3 components.
    """
    def __init__(self, z=[0., 0., 0.], curvature=1, Ra=50.):
        self.__curvature = float(curvature)
        OpticalElement.__init__(self, z, Ra)
    
    def centre(self):
        """Return the centre of curvature of the surface as a 3-dimensional
        numpy array. If curvature is zero, return the axial intersection.
        """
        if self.__curvature == 0:
                return self._z                                               
        Rc = (self.__curvature)**(-1)
        u = np.array([0., 0., 1.])
        return self._z + Rc*u
        
    def intercept(self, ray):
        """Return the point on the surface as a 3-dimensional numpy array where 
        a given ray is intercepted by the SphericalReflection object.
        
        Parameters: ray -- Ray object.
        """
        length = 0
        if ray._terminate is True:
            return None
        if self.__curvature == 0:                                                
            shift_p = ray.p() - np.array([ray.p()[0], 0., 0.])\
                      - np.array([0., ray.p()[1], 0.])
            theta = np.arccos(np.dot(np.array([0, 0, 1]), unit(ray.k())))
            r = shift_p - self._z
            length = abs(la.norm(r)*(np.cos(theta))**(-1))
            intercept = ray.p() + length*unit(ray.k())
            if np.sqrt(intercept[0]**2. + intercept[1]**2.) > self._Ra:
                return None 
            else:
                return intercept
        else: 
            Rc = self.__curvature**(-1)
            r = ray.p() - self.centre()
        r_dot_k = np.dot(r, unit(ray.k()))
        if (r_dot_k)**2. - np.dot(r, r) + Rc**2. < 0:                             
            return None                                                         
        length_1 = -r_dot_k + ((r_dot_k)**2. - np.dot(r, r) + Rc**2.)**(0.5)
        length_2 = -r_dot_k - ((r_dot_k)**2. - np.dot(r, r) + Rc**2.)**(0.5)
        if self.__curvature > 0:                                                 
            length = min(abs(length_1), abs(length_2))
        elif self.__curvature < 0:                                               
            if la.norm(ray.p() - self.centre()) < abs(Rc):                          
                if np.dot(ray.p() - self.centre(),
                          np.array([0., 0., 1.])) >= 0:
                    length = min(abs(length_1), abs(length_2))
                else:
                    length = max(abs(length_1), abs(length_2))
            else:
                length = max(abs(length_1), abs(length_2))
        intercept = ray.p() + length*unit(ray.k())
        if np.sqrt(intercept[0]**2 + intercept[1]**2) > self._Ra:
            return None                                                         
        else:                                                                   
            return intercept
            
    def propagate_ray(self, ray):
        """Propagate a given ray from its current point and reflect the ray
        in the direction calculated by the reflect function. Change the ray's 
        current direction and append intercept point to ray's internal list
        of points.
        
        Parameters: ray -- Ray object.
        """
        if self.intercept(ray) is None:                                         
            ray.terminate()
            print "Ray did not propagate."
            return None
        if self.__curvature == 0.0:
            surface_normal = -unit(self._z)
        else: 
            surface_normal = unit(self.intercept(ray) - self.centre())
        new_direction = reflect(ray.k(), surface_normal)
        ray.append(self.intercept(ray), new_direction)
        