# -*- coding: utf-8 -*-

###############################################################################
###############################################################################
##                                                                           ##
##    FILE DESCRIPTION:                                                      ##
##                                                                           ##
##    Function to solve for the orbit position, velocity and true anomaly.   ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 02-May-2021 00:53 AM (+8 GMT)                            ##
##    Last modified 20-May-2021 10:14 AM (+8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import math
import numpy as np

from source import anomaly
from source import rotation


class Spacecraft:
    
    def __init__(self, elements=None, posvel=None, mass=1, area=1, Cd=2.2):
        '''Initialise a spacecraft object defined by its 6 states, mass, drag
        area, and coefficient of drag.

        Parameters
        ----------
        elements : list, optional
            List of 6 Keplerian elements (a, e, i, R, w, M). Declare
            only either the elements or the position and velocity vectors.
        posvel : list, optional
            List of 6 Cartesian states (px, py, pz, vx, vy, vz). Declare 
            only either the elements or the position and velocity vectors.
        mass : float, optional
            Drag mass of spacecraft (kg). The default is 0.
        area : float, optional
            Drag area of spacecraft (m^2). The default is 0.
        Cd : float, optional
            Dimensionless drag coefficient. The default is 2.2.
            
        Attributes
        ----------
        
        '''
        
        # G * Earth Mass (km**3/s**2)
        GM = 398600.4418
        
        self.Ms = mass # kg
        self.Ar = area # m^2
        self.Cd = Cd   # unitless
        
        # Revert to defaults if no constructor arguments are specified.
        if elements is None and posvel is None:
            
            self.px = 0.0 # Inertial position in X-axis (km)
            self.py = 0.0 # Inertial position in Y-axis (km)
            self.pz = 0.0 # Inertial position in Z-axis (km)
            self.vx = 0.0 # Inertial velocity in X-axis (km/s)
            self.vy = 0.0 # Inertial velocity in Y-axis (km/s)
            self.vz = 0.0 # Inertial velocity in Z-axis (km/s)
            
            self.a  = 0.0 # Semi-major axis (m)
            self.e  = 0.0 # Eccentricity (unitless)
            self.i  = 0.0 # Inclination (radians)
            self.w  = 0.0 # Arg of Periapsis (radians)
            self.R  = 0.0 # Right Ascension (radians)
            self.M  = 0.0 # Mean Anomaly (radians)
            self.U  = 0.0 # Mean Arg Latitude (radians)
            self.nu = 0.0 # True anomaly (radians)
            self.n  = 0.0 # Mean motion (radians/second)
        
        # Compute the positions and velocities if elements are specified.
        elif elements is not None:
            
            try:
                # Get the orbital elements.
                self.a = elements[0] # Semi-major axis (m)
                self.e = elements[1] # Eccentricity (unitless)
                self.i = elements[2] # Inclination (radians)
                self.w = elements[3] # Arg of Periapsis (radians)
                self.R = elements[4] # Right Ascension (radians)
                self.M = elements[5] # Mean Anomaly (radians)
                
                # Compute the positions and velocities.
                px, py, pz, vx, vy, vz, nu = self.update_by_elements(self.a,
                                                                     self.e,
                                                                     self.i,
                                                                     self.w,
                                                                     self.R,
                                                                     self.M)
                
                # Store the positions and velocities
                self.px = px # Inertial position in X-axis (km)
                self.py = py # Inertial position in Y-axis (km)
                self.pz = pz # Inertial position in Z-axis (km)
                self.vx = vx # Inertial velocity in X-axis (km/s)
                self.vy = vy # Inertial velocity in Y-axis (km/s)
                self.vz = vz # Inertial velocity in Z-axis (km/s)
                self.nu = nu     # True anomaly (radians)
                
                # Update mean argument of latitude
                self.U = self.M + self.w
                
                # Update the mean motion (radians/second).
                self.n = math.sqrt( GM / (self.a**3) )
            
            except TypeError:
                print("TypeError: Elements must be a list [a,e,i,w,R,M]! \n")
            except IndexError:
                print("IndexError: Elements list length must be 6! \n")
            except:
                print("Unknown error occurred. Printing constructor args:")
                print(elements, posvel, mass, area, Cd, " \n")
        
        # Update elements if the positions and velocities are specified.
        elif posvel is not None:
            
            try:
                # Get the positions and velocities
                self.px = posvel[0] # Inertial position in X-axis (km)
                self.py = posvel[1] # Inertial position in Y-axis (km)
                self.pz = posvel[2] # Inertial position in Z-axis (km)
                self.vx = posvel[3] # Inertial velocity in X-axis (km/s)
                self.vy = posvel[4] # Inertial velocity in Y-axis (km/s)
                self.vz = posvel[5] # Inertial velocity in Z-axis (km/s)
                
                # Compute the orbital elements.
                a, e, i, w, R, M, nu = self.update_by_posvel(self.px,
                                                             self.py,
                                                             self.pz,
                                                             self.vx,
                                                             self.vy,
                                                             self.vz)
                
                self.a  = a  # Semi-major axis (m)
                self.e  = e  # Eccentricity (unitless)
                self.i  = i  # Inclination (radians)
                self.w  = w  # Arg of Periapsis (radians)
                self.R  = R  # Right Ascension (radians)
                self.M  = M  # Mean Anomaly (radians)
                self.nu = nu # True anomaly (radians)
                
                # Update mean argument of latitude
                self.U = self.M + self.w
                
                # Update the mean motion.
                self.n = math.sqrt( GM / (self.a**3) )
            
            except TypeError:
                print("TypeError: Pos-Vel must be a list of 6 states!")
            except IndexError:
                print("IndexError: Pos-Vel list length must be 6!")
            except:
                print("Unknown error occurred. Printing constructor args:")
                print(elements, posvel, mass, area, Cd)
        
    def update_by_elements(self, a, e, i, w, R, M):
        '''Method that returns an inertial position vector (1x3 NumPy 
        array), an inertial velocity vector (1x3), and a true anomaly value
        (float), when ingesting six osculating Keplerian orbit elements.
        Reference: "Satellite Orbits" by Oliver Montenbruck, Chapter 2.2.3.
        
        Parameters
        ----------
        a : float
            Semi-major axis (km)
        e : float
            Eccentricity (unit-less)
        i : float
            Inclination (rad)
        w : float
            Argument of Perigee (rad)
        R : float
            Right Angle of Asc Node (rad)
        M : float
            Mean Anomaly (rad)
    
        Returns
        -------
        pos : numpy.ndarray
            Inertial position vector (1x3 vector, km)
        vel : numpy.ndarray
            Inertial velocity vector (1x3 vector, km/s)
        nu  : numpy.float64
            True anomaly (float, rad)
        
        '''
        
        # G * Earth Mass (km**3/s**2)
        GM = 398600.4418
        
        # The general flow of the program, is to first solve for the radial
        # position and velocity (in the inertial frame) via Kepler's equation.
        # Thereafter, we obtain the inertial coordinates in the Hill frame,
        # by performing a 3-1-3 Euler Angle rotation using an appropriate DCM.
        
        # First, let us solve for the eccentric anomaly.
        eccAnom = anomaly.M2E(M,e)
        
        # With the eccentric anomaly, we can solve for position and velocity
        # in the local orbital frame, using the polar equation for an ellipse.
        pos_X = a * ( np.cos(eccAnom) - e)
        pos_Y = a * np.sqrt( 1 - e**2 ) * np.sin(eccAnom)
        pos_norm = np.sqrt( pos_X**2 + pos_Y**2 )
        vel_const = np.sqrt( GM * a ) / pos_norm
        vel_X = vel_const * ( -1 * np.sin(eccAnom) )
        vel_Y = vel_const * ( np.sqrt( 1 - e**2 ) * np.cos(eccAnom) )
        
        # To perform the conversion from local orbit plane to an ECI frame, we
        # need perform the 313 Euler angle rotation in the following sequence:
        # Right Angle of Ascending Node -> Inclination -> Argument of Latitude
        # Now, let us get us the DCM that converts to the hill-frame.
        
        DCM_HN = rotation.dcmZ(w) @ rotation.dcmX(i) @ rotation.dcmZ(R)
        
        # Notice that the hill frame computation does not include a rotation
        # of the true anomaly, and that's because the true anomaly has already
        # been accounted for when computing pos_X and pos_Y using information 
        # from the eccentric anomaly. Including true anomaly in the DCM 
        # rotation would double-count that anomaly rotation.
        
        # The current coordinates are in the local hill frame, and thus 
        # conversion from hill to inertial would be the transpose of HN.
        DCM_NH = np.transpose(DCM_HN)
        
        # With the hill frame, we can now convert it to the ECI frame.
        pos = DCM_NH @ np.array([ pos_X, pos_Y, 0.0 ]).T
        vel = DCM_NH @ np.array([ vel_X, vel_Y, 0.0 ]).T
        
        # Finally, let us not forget to compute the true anomaly.
        nu = np.arctan2( pos_Y, pos_X )
        
        # Update self values.
        self.px, self.py, self.pz = pos[0], pos[1], pos[2]
        self.vx, self.vy, self.vz = vel[0], vel[1], vel[2]
        self.a, self.e, self.i, self.w, self.R, self.M  = a, e, i, w, R, M
        self.nu = nu
        
        # Update the mean motion.
        self.n = math.sqrt( GM / (self.a**3) )
        
        # Position 1x3 (km), velocity 1x3 (km/s), true anomaly (rad)
        return self.px, self.py, self.pz, self.vx, self.vy, self.vz, self.nu
    
    def update_by_posvel(self, px, py, pz, vx, vy, vz):
        '''Method that returns six osculating Keplerian orbit elements
        (a, e, i, w, R, M), and ingests in a 1x3 position vector (inertial)
        and a 1x3 velocity vector (inertial). Reference: "Satellite Orbits"
        by Oliver Montenbruck, Chapter 2.2.4.
        
        Parameters
        ----------
        pos : numpy.ndarray
            Inertial position vector (1x3 vector, km)
        vel : numpy.ndarray
            Inertial velocity vector (1x3 vector, km/s)
        
        Returns
        -------
        a : float
            Semi-major axis (km)
        e : float
            Eccentricity (unit-less)
        i : float
            Inclination (rad)
        w : float
            Argument of Perigee (rad)
        R : float
            Right Angle of Asc Node (rad)
        M : float
            Mean Anomaly (rad)
        nu : float
            True Anomaly (rad)
        
        '''
        
        GM  = 398600.4418
        pos = np.array([px, py, pz])
        vel = np.array([vx, vy, vz])
        
        # First, compute the semi-major axis.
        r = np.linalg.norm(pos)
        a = 1 / ( (2/r) - ( ( (np.linalg.norm(vel))**2 ) / GM ) )
        
        # Second, compute the angular momentum vector of the orbit.
        H      = np.cross(pos,vel)
        H_norm = np.linalg.norm(H)
        H_hat  = H / H_norm
        Wx     = H_hat[0]
        Wy     = H_hat[1]
        Wz     = H_hat[2]
        
        # Third, from the normalised angular momentum, derive the inclination.
        i = np.arctan2( math.sqrt( Wx**2 + Wy**2 ), Wz )
        
        # Fourth, from the normalised angular momentum, derive the RAAN.
        R = np.arctan2( Wx, -1*Wy )
        
        # Fifth, compute the semi-latus rectum.
        p = ( H_norm**2 / GM )
        
        # Sixth, fetch the mean motion.
        n = self.n
        
        # Seventh, assuming an elliptical orbit, compute the eccentricity.
        e = math.sqrt( 1 - (p/a) )
        
        # Eighth, compute the eccentric anomaly.
        E = np.arctan2( ( (np.dot(pos,vel) / ((a**2)*n)) ), (1 - r/a) )
        
        # Ninth, we can compute the mean anomaly using Kepler's equation.
        M = E - e*math.sin(E)
        
        # Tenth, the argument of latitude is computed.
        U = np.arctan2( pos[2], ( pos[1]*Wx - pos[0]*Wy ) )
        
        # Eleventh, the true anomaly is computed.
        nu = np.arctan2( math.sin(E)*math.sqrt(1-e**2), math.cos(E)-e )
        
        # Twelfth, the argument of perigee is computed.
        w = U - nu
        
        # Update self values.
        self.px, self.py, self.pz = px, py, pz
        self.vx, self.vy, self.vz = vx, vy, vz
        self.a, self.e, self.i, self.w, self.R, self.M  = a, e, i, w, R, M
        self.nu = nu
        
        # Update mean argument of latitude
        self.U = M + w
        
        # Update the mean motion.
        self.n = math.sqrt( GM / (self.a**3) )
        
        return self.a, self.e, self.i, self.w, self.R, self.M, self.nu
    
    def twobody_propagate(self, tstep):
        
        # Wrap the mean anomaly about +/- 180 degrees.
        self.M = (self.M + (self.n * tstep) + math.pi) % (2*math.pi) - math.pi
        self.update_by_elements(self.a, self.e, self.i, self.w, self.R, self.M)
        