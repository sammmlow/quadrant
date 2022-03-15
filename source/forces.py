# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 22:08:12 2021

@author: sammm
"""

import numpy as np
from source import atmosphere

def forces( pos, vel, sc ):
    '''Computation of the total inertial acceleration as a 1x3 vector, from
    Earth's gravity; optionally the J2 perturbation force, and drag force via
    the US Standard Atmosphere 1976.
    
    Parameters
    ----------
    sc : numpy.ndarray
        xxx
        
    Returns
    -------
    acceleration : numpy.ndarray
        Inertial frame acceleration vector (1x3) of the spacecraft (km/s^2)
    
    '''
    
    # Retrieve all parameters from the spacecraft.
    GM = sc.GM
    Cd = sc.Cd
    
    # Define all constants
    RE = 6378.140     # Earth equatorial radius (km)
    GM = 398600.4418  # G * Earth Mass (km**3/s**2)
    J2 = 1.0826267e-3 # J2 constant
    
    # Get the radial distance of the satellite.
    R = np.linalg.norm( pos ) # km
    V = np.linalg.norm( vel ) # km/s
    
    # Initialise the acceleration vector.
    acceleration = np.zeros(3)
    
    # Compute the two-body gravitational force by Earth.
    if sc.forces['Earth Twobody'] == True:
        acceleration += ( -1 * GM * pos ) / ( R**3 )
    
    # Include the additional J2 acceleration vector if necessary.
    if sc.forces['Earth Oblate J2'] == True:
        R_J2 = 1.5 * J2 * GM * ((RE**2)/(R**5))
        zRatio = (pos[2]/R)**2
        oblate_x = R_J2 * pos[0] * (5 * zRatio-1)
        oblate_y = R_J2 * pos[1] * (5 * zRatio-1)
        oblate_z = R_J2 * pos[2] * (5 * zRatio-3)
        acceleration += np.array([oblate_x, oblate_y, oblate_z])
    
    # Include the additional drag acceleration if necessary.
    if sc.forces['Earth Atmos Drag'] == True:
        areaMassRatio = sc.area / sc.mass # m**2/kg
        dragDensity = atmosphere.density( (R - RE) ) # kg/m**3
        dragAccel = 0.5 * Cd * dragDensity * areaMassRatio * ((V*1000)**2)
        acceleration -= dragAccel * ( vel / V ) / 1000
    
    # Acceleration vector is in km/s**2
    return acceleration
