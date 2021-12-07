# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 22:08:12 2021

@author: sammm
"""

import numpy as np
from source import atmos

def acceleration( pos, vel, Cd, Ar, Ms, fJ, fD ):
    '''Computation of the total acceleration as a 1x3 vector, primarily due to
    Earth gravity, and optionally (if fJ == 1) the J2 perturbation force, and
    optionally (if fD == 1) the drag force via US Standard Atmosphere 1976.
    
    Parameters
    ----------
    pos : numpy.ndarray
        Inertial frame position vector (1x3) of the spacecraft (km)
    vel : numpy.ndarray
        Inertial frame velocity vector (1x3) of the spacecraft (km/s)
    Cd : float
        Drag coefficient of the spacecraft
    Ar : float
        Drag area of the spacecraft (m^2)
    Ms : float
        Mass of the spacecraft (kg)
    fJ : bool
        Flag to toggle J2 perturbation (True to toggle on)
    fD : bool
        Flag to toggle atmospheric drag (True to toggle on)
    
    Returns
    -------
    acceleration : numpy.ndarray
        Inertial frame acceleration vector (1x3) of the spacecraft (km/s^2)
    
    '''
    
    # Define all constants
    RE = 6378.140     # Earth equatorial radius (km)
    GM = 398600.4418  # G * Earth Mass (km**3/s**2)
    J2 = 1.0826267e-3 # J2 constant
    
    # Get the radial distance of the satellite.
    R = np.linalg.norm( pos ) # km
    V = np.linalg.norm( vel ) # km/s
    
    # Compute the two-body gravitational force by Earth.
    acceleration = ( -1 * GM * pos ) / ( R**3 )
    
    # Include the additional J2 acceleration vector if necessary.
    if fJ == True:
        R_J2 = 1.5 * J2 * GM * ((RE**2)/(R**5))
        zRatio = (pos[2]/R)**2
        oblate_x = R_J2 * pos[0] * (5 * zRatio-1)
        oblate_y = R_J2 * pos[1] * (5 * zRatio-1)
        oblate_z = R_J2 * pos[2] * (5 * zRatio-3)
        acceleration = acceleration + np.array([oblate_x, oblate_y, oblate_z])
    
    # Include the additional drag acceleration if necessary.
    if fD == True:
        areaMassRatio = Ar / Ms; # m**2/kg
        dragDensity = atmos.density( (R - RE) ) # kg/m**3
        dragAccel = 0.5 * Cd * dragDensity * areaMassRatio * ( (V*1000)**2 )
        acceleration = acceleration - ( dragAccel * ( vel / V ) / 1000 )
    
    # Acceleration vector is in km/s**2
    return acceleration
