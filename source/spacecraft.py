# -*- coding: utf-8 -*-

###############################################################################
###############################################################################
##                                                                           ##
##      ___  _   _   __   ____  ____   __   _   _ _____                      ##
##     / _ \| | | | /  \ |  _ \| __ \ /  \ | \ | |_   _|                     ##
##    ( |_| ) |_| |/ /\ \| |_| | -/ // /\ \|  \| | | |                       ##
##     \_  /|_____| /--\ |____/|_|\_\ /--\ |_\___| |_|                       ##
##       \/                                               v 0.0              ##
##                                                                           ##
##    Single spacecraft class.                                               ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 02-May-2021 00:53 AM (+8 GMT)                            ##
##    Last modified 15-Mar-2022 13:05 AM (-8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import math
import random
import datetime
import numpy as np

from source import anomaly
from source import rotation
from source import integrate
from source.attitudes import QTR, CRP, MRP

PI = 3.141592653589793
R2D = 180.0/PI # Radians to degrees
D2R = PI/180.0 # Degrees to radians

# Random names for spacecraft from (thanks FTL)
names = ['Arries', 'Elnubnub', 'Ken', 'Randy', 'Kirkner', 'Lombard', 'Feras',
         'Yeoz', 'Ryan', 'GM Faux', 'O Williams', 'Yevon Si', 'Maloney',
         'Will', 'Maxwell', 'Thomas', 'Bloch', 'Rynhart', 'AJ Hager',
         'Stephen', 'Eckman', 'Jack', 'Shirai', 'Julian', 'Rafn', 'Markand',
         'Nelis', 'Chris', 'Malott', 'Davion', 'Caldwell', 'Kusy', 'Lagardi',
         'Vincent', 'Kevin', 'Fish', 'Turing', 'Stelly', 'Teldarin']

class Spacecraft():
    '''
    Initialise a spacecraft object with state attributes, attitude attributes,
    angular velocities, and area-mass parameters.
    
    Attributes
    ----------
    
    elements : list (optional)
        List of 6 Keplerian elements [a, e, i, R, w, M].
        
    states : list (optional)
        List of 6 Cartesian states [px, py, pz, vx, vy, vz].
        
    mass : float, optional
        Drag mass of spacecraft (kg). The default is 1.
        
    area : float, optional
        Drag area of spacecraft (m^2). The default is 0.
        
    Cd : float, optional
        Dimensionless drag coefficient. The default is 2.2.
        
    name : str, optional
        Custom name for the spacecraft object.
        
    GM : float, optional
        Product of planet mass and gravitational constant (km**3/s**2)
        
    '''
    
    # Set default parameters upon initialization, in this order.
    def __init__( self, elements=None, states=None, epoch=None, forces=None, 
                  torque=None, attBN=None, ohmBN=None, name=None, mass=None,
                  inertia=None, area=None, Cd=None, GM=None ):
                 
                 
        
        # State vector units are:
        # elements: [ km, none, deg, deg, deg, deg ]
        # states:   [ km, km, km, km/s, km/s, km/s ]
        
        # Initialize non-state parameters.
        self.epoch = epoch
        self.forces = forces
        self.torque = torque
        self.attBN = attBN
        self.ohmBN = ohmBN
        self.attBR = attBN
        self.ohmBR = ohmBN
        self.name = name
        self.mass = mass
        self.inertia = inertia
        self.area = area
        self.Cd = Cd
        self.GM = GM
        self.attIntgErr = None
        
        # If the elements are defined, then initialize both elements and the
        # cartesian states from `elements`. Else, initialize it via `states`.
        if elements is not None and states is None:
            self._resetstates()
            self.elements = elements
        elif elements is None and states is not None:
            self._resetstates()
            self.states = states
        else:
            self._resetstates()
            self.elements = elements
            self.states = states
    
    # Reset spacecraft state attributes.
    def _resetstates(self):
        self.__dict__['a']  = 0.0
        self.__dict__['e']  = 0.0
        self.__dict__['i']  = 0.0
        self.__dict__['w']  = 0.0
        self.__dict__['R']  = 0.0
        self.__dict__['M']  = 0.0
        self.__dict__['px'] = 0.0
        self.__dict__['py'] = 0.0
        self.__dict__['pz'] = 0.0
        self.__dict__['vx'] = 0.0
        self.__dict__['vy'] = 0.0
        self.__dict__['vz'] = 0.0
        self.__dict__['nu'] = 0.0
        self.__dict__['mu'] = 0.0
        self.__dict__['n']  = 0.0
        self.__dict__['T']  = 0.0
        
    # Constrain the attribute values and types strictly in this order.
    def __setattr__(self, key, value):
        
        # __setattr__ for 'epoch' attribute. Datetime object.
        if key == 'epoch':
            if value is None:
                now = datetime.datetime.now().replace(microsecond=0)
                self.__dict__[key] = now
            else:
                if type(value) == datetime.datetime:
                    self.__dict__[key] = value
                else:
                    raise ValueError("Epoch must be a datetime object!")
        
        # __setattr__ for 'forces' attribute. Dictionary of forces toggled.
        if key == 'forces':
            if value is None:
                self.__dict__[key] = {'Earth Twobody':True,
                                      'Earth Oblate J2':False,
                                      'Earth Atmos Drag':False,
                                      'Earth Solar Rad':False,
                                      'Moon Third-Body':False,
                                      'Sun Third-Body':False}
            else:
                if type(value) is dict:
                    print('Toggling the current forces:')
                    print('----------------------------')
                    for element in value:
                        print(element, 'toggled to', str(value[element]))
                    self.__dict__[key] = value
                else:
                    print('Error in force model, printing:')
                    print(value)
                    raise TypeError("Force models must be in a dictionary!")
        
        # __setattr__ for 'torque' attribute. Dictionary of forces toggled.
        if key == 'torque':
            if value is None:
                self.__dict__[key] = np.array([0.0,0.0,0.0])
            else:
                if len(value) == 3:
                    self.__dict__[key] = np.array(value)
                else:
                    self.__dict__[key] = np.array([0.0,0.0,0.0])
                    raise ValueError("Torque length not 3!")
        
        # __setattr__ for 'attBN' attribute. Body to inertial attitude.
        if key == 'attBN':
            if value is None:
                self.__dict__[key] = QTR()
            else:
                if len(value) == 4 and value.strID() == 'QTR':
                    self.__dict__[key] = value
                elif len(value) == 3:
                    if value.strID() == 'CRP':
                        self.__dict__[key] = value
                    elif value.strID() == 'MRP':
                        self.__dict__[key] = value
                    else:
                        raise TypeError("Attitude BN type unknown!")
                else:
                    raise ValueError("Attitude BN vector length incorrect!")
        
        # __setattr__ for 'ohmBN' attribute.
        if key == 'ohmBN':
            if value is None:
                self.__dict__[key] = np.zeros(3)
            else:
                if len(value) == 3:
                    self.__dict__[key] = np.array( value )
                else:
                    raise ValueError("Omega BN vector length incorrect!")
        
        # __setattr__ for 'attBR' attribute.
        if key == 'attBR':
            if value is None:
                self.__dict__[key] = QTR()
            else:
                if len(value) == 4 and value.strID() == 'QTR':
                    self.__dict__[key] = value
                elif len(value) == 3:
                    if value.strID() == 'CRP':
                        self.__dict__[key] = value
                    elif value.strID() == 'MRP':
                        self.__dict__[key] = value
                    else:
                        raise TypeError("Attitude BR type unknown!")
                else:
                    raise ValueError("Attitude BR length incorrect!")
        
        # __setattr__ for 'ohmBR' attribute.
        if key == 'ohmBR':
            if value is None:
                self.__dict__[key] = np.zeros(3)
            else:
                if len(value) == 3:
                    self.__dict__[key] = np.array( value )
                else:
                    raise ValueError("Omega BR vector length incorrect!")
        
        # __setattr__ for 'attIntgErr' attribute.
        if key == 'attIntgErr':
            if value is None:
                self.__dict__[key] = np.zeros(3)
            else:
                if len(value) == 3:
                    self.__dict__[key] = value
                else:
                    self.__dict__[key] = np.zeros(3)
                    raise ValueError("Attitude integral error not length 3!")
        
        # __setattr__ for 'name' attribute.
        if key == 'name':
            if value is None:
                self.__dict__[key] = random.choice(names)
            else:
                self.__dict__[key] = str(value)
        
        # __setattr__ for 'mass' attribute.
        if key == 'mass':
            if value is None:
                self.__dict__[key] = 100.0
            else:
                self.__dict__[key] = float(value)
        
        # __setattr__ for 'mass' attribute.
        if key == 'inertia':
            if value is None:
                self.__dict__[key] = np.diag([10,10,10])
            else:
                if np.shape(value) == (3,3):
                    self.__dict__[key] = value
                    if sum(sum(value)) != np.trace(value):
                        print("Warning, inertia tensor is not diagonal!")
                else:
                    raise ValueError("Inertia tensor must be 3x3 diagonal!")
        
        # __setattr__ for 'area' attribute.
        if key == 'area':
            if value is None:
                self.__dict__[key] = 1.0
            else:
                self.__dict__[key] = value
        
        # __setattr__ for 'Cd' attribute.
        if key == 'Cd':
            if value is None:
                self.__dict__[key] = 2.2
            else:
                self.__dict__[key] = value
        
        # __setattr__ for 'GM' attribute. Product of the universal gravity
        # constant and the central body's mass. Defaults to Earth if None.
        # Units must be in (km**3/s**2)
        if key == 'GM':
            if value is None:
                self.__dict__[key] = 398600.4418 # Earth default (km**3/s**2)
            else:
                self.__dict__[key] = value
                
        # __setattr__ for 'elements' attribute. Array of 6 Keplerian elements.
        if key == 'elements':
            if value is None:
                self.__dict__[key] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                self._resetstates()
            else:
                if len(value) == 6:
                    a = value[0]
                    e = value[1]
                    i = (D2R * value[2]) % (PI)
                    w = ((((D2R * value[3]) + PI) % (2*PI)) - PI)
                    R = ((((D2R * value[4]) + PI) % (2*PI)) - PI)
                    M = ((((D2R * value[5]) + PI) % (2*PI)) - PI)
                    self.__dict__[key] = [a,e,i,w,R,M]
                    self.__dict__['a'] = a
                    self.__dict__['e'] = e
                    self.__dict__['i'] = i
                    self.__dict__['w'] = w
                    self.__dict__['R'] = R
                    self.__dict__['M'] = M
                    self._update_states()
                else:
                    raise ValueError("Wrong length of orbit elements!")
        
        # __setattr__ for 'states' attribute. Units in km and km/s.
        if key == 'states':
            if value is None:
                self.__dict__[key] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                self._resetstates()
            else:
                if len(value) == 6:
                    self.__dict__[key] = value
                    self.__dict__['px'] = value[0]
                    self.__dict__['py'] = value[1]
                    self.__dict__['pz'] = value[2]
                    self.__dict__['vx'] = value[3]
                    self.__dict__['vy'] = value[4]
                    self.__dict__['vz'] = value[5]
                    self._update_elements()
                else:
                    raise ValueError("Wrong length of spacecraft states!")
        
        # __setattr__ for 'px' attribute. Inertial position X, units in km.
        if key == 'px':
            self.__dict__[key] = value
            self.__dict__['states'][0] = value
            self._update_elements()
        
        # __setattr__ for 'py' attribute. Inertial position Y, units in km.
        if key == 'py':
            self.__dict__[key] = value
            self.__dict__['states'][1] = value
            self._update_elements()
        
        # __setattr__ for 'pz' attribute. Inertial position Z, units in km.
        if key == 'pz':
            self.__dict__[key] = value
            self.__dict__['states'][2] = value
            self._update_elements()
        
        # __setattr__ for 'vx' attribute. Inertial velocity X, units in km/s.
        if key == 'vx':
            self.__dict__[key] = value
            self.__dict__['states'][3] = value
            self._update_elements()
        
        # __setattr__ for 'vy' attribute. Inertial velocity Y, units in km/s.
        if key == 'vy':
            self.__dict__[key] = value
            self.__dict__['states'][4] = value
            self._update_elements()
        
        # __setattr__ for 'vz' attribute. Inertial velocity Z, units in km/s.
        if key == 'vz':
            self.__dict__[key] = value
            self.__dict__['states'][5] = value
            self._update_elements()
        
        # __setattr__ for 'a' attribute. Element semi-major axis (km).
        if key == 'a':
            self.__dict__[key] = value
            self.__dict__['elements'][0] = value
            self._update_states()
            self.__dict__['n'] = (self.GM / (value**3))**0.5
            self.__dict__['T'] = 2 * PI * ((value**3)/self.GM)**0.5
        
        # __setattr__ for 'e' attribute. Element eccentricity (unit-less).
        if key == 'e':
            self.__dict__[key] = value
            self.__dict__['elements'][1] = value
            self._update_states()
            if value >= 1.0:
                raise ValueError("Eccentricity must be less than 1!")
        
        # __setattr__ for 'i' attribute. Element inclination (rad).
        if key == 'i':
            self.__dict__[key] = (value * D2R) % PI
            self.__dict__['elements'][2] = (value * D2R)
            self._update_states()
        
        # __setattr__ for 'w' attribute. Element argument of periapsis (rad).
        if key == 'w':
            self.__dict__[key] = (((value * D2R) + PI) % (2*PI)) - PI
            self.__dict__['elements'][3] = (((value * D2R) + PI) % (2*PI)) - PI
            self._update_states()
        
        # __setattr__ for 'R' attribute. Element right ascending node (rad)
        if key == 'R':
            self.__dict__[key] = (((value * D2R) + PI) % (2*PI)) - PI
            self.__dict__['elements'][4] = (((value * D2R) + PI) % (2*PI)) - PI
            self._update_states()
        
        # __setattr__ for 'M' attribute. Element mean anomaly (rad).
        if key == 'M':
            self.__dict__[key] = (((value * D2R) + PI) % (2*PI)) - PI
            self.__dict__['elements'][5] = (((value * D2R) + PI) % (2*PI)) - PI
            self._update_states()
        
        # __setattr__ for 'nu' attribute. True anomaly (rad).
        if key == 'nu':
            self.__dict__[key] = (((value * D2R) + PI) % (2*PI)) - PI
            self.M = R2D * (((anomaly.V2M((value*D2R),self.e)+PI)%(2*PI))-PI)
        
        # __setattr__ for 'mu' attribute. Mean argument of latitude (rad).
        if key == 'mu':
            self.__dict__[key] = ((value + PI) % (2*PI)) - PI
            self.nu = R2D * ((((value * D2R) - self.w + PI) % (2*PI)) - PI)
        
        # __setattr__ for 'n' attribute. Mean motion (rad/s). Constrained by
        # semi-major axis value (km).
        if key == 'n':
            self.__dict__[key] = value
            self.__dict__['a'] = (self.GM / value**2)**(1/3)
            self.__dict__['T'] = 2 * PI / value
        
        # __setattr__ for 'T' attribute. Keplerian period (s). Constrained by
        # semi-major axis value (km).
        if key == 'T':
            self.__dict__[key] = value
            self.__dict__['n'] = 2 * PI / value
            self.__dict__['a'] = (self.GM * (value / (2*PI))**2 )**(1/3)
    
    def __repr__(self):
        return 'Spacecraft ' + str(self.name)
    
    def status(self):
        
        # Print out basic spacecraft information.
        out  = 'BASIC INFORMATION \n'
        out += '----------------- \n'
        out += 'Current Date-Time   : ' + str(self.epoch) + ' \n'
        out += 'Spacecraft Name     : ' + str(self.name)  + ' \n'
        out += 'Spacecraft Wet Mass : ' + str(self.mass)  + ' kg\n'
        out += 'Drag Area Effective : ' + str(self.area)  + ' m**2\n'
        out += 'Drag Coefficient    : ' + str(self.Cd)    + ' \n'
        out += 'Gravitational Const : ' + str(self.GM)    + ' km**3/s**2\n'
        out += '\n'
        
        # Print out inertial Cartesian states.
        out += 'INERTIAL STATES \n'
        out += '--------------- \n'
        out += 'Inertial Position X : '+'{:.6f}'.format(self.px)    +' km\n'
        out += 'Inertial Position Y : '+'{:.6f}'.format(self.py)    +' km\n'
        out += 'Inertial Position Z : '+'{:.6f}'.format(self.pz)    +' km\n'
        out += 'Inertial Velocity X : '+'{:.6f}'.format(self.vx)    +' km/s\n'
        out += 'Inertial Velocity Y : '+'{:.6f}'.format(self.vy)    +' km/s\n'
        out += 'Inertial Velocity Z : '+'{:.6f}'.format(self.vz)    +' km/s\n'
        out += '\n'
        
        # Print out osculating orbital elements and orbital parameters.
        out += 'OSCULATING ORBIT ELEMENTS AND PARAMETERS \n'
        out += '---------------------------------------- \n'
        out += 'Osc Semi-Major Axis : '+'{:.6f}'.format(self.a)     +' km\n'
        out += 'Osc Eccentricity    : '+'{:.6f}'.format(self.e)     +' \n'
        out += 'Osc Inclination     : '+'{:.6f}'.format(self.i*R2D) +' deg\n'
        out += 'Osc Arg. Periapsis  : '+'{:.6f}'.format(self.w*R2D) +' deg\n'
        out += 'Osc Right Asc Node  : '+'{:.6f}'.format(self.R*R2D) +' deg\n'
        out += 'Osc Mean Anomaly    : '+'{:.6f}'.format(self.M*R2D) +' deg\n'
        out += 'Osc True Anomaly    : '+'{:.6f}'.format(self.nu*R2D)+' deg\n'
        out += 'Osc Arg. Latitude   : '+'{:.6f}'.format(self.mu*R2D)+' deg\n'
        out += 'Mean Motion         : '+'{:.6f}'.format(self.n)     +' rad/s\n'
        out += 'Keplerian Period    : '+'{:.6f}'.format(self.T)     +' s\n'
        out += '\n'
        
        # Print out attitude coordinates
        out += 'ATTITUDE COORDINATES AND ANGULAR RATES \n'
        out += '-------------------------------------- \n'
        if self.attBN.strID() == 'QTR':
            out += 'Attitude Type BN    : Quaternion \n'
        elif self.attBN.strID() == 'CRP':
            out += 'Attitude Type BN    : Classical Rodrigues Parameters \n'
        elif self.attBN.strID() == 'MRP':
            out += 'Attitude Type BN    : Modified Rodrigues Parameters \n'
        else:
            out += 'Attitude Type BN    : ERRONEOUS COORDINATE \n'
        out += 'Attitude Coord BN   : '+str(np.around(self.attBN,5))+' \n'
        out += 'Angular Velocity BN : '+str(np.around(self.ohmBN,5))+' rad/s\n'
        if self.attBR.strID() == 'QTR':
            out += 'Attitude Type BR    : Quaternion \n'
        elif self.attBR.strID() == 'CRP':
            out += 'Attitude Type BR    : Classical Rodrigues Parameters \n'
        elif self.attBR.strID() == 'MRP':
            out += 'Attitude Type BR    : Modified Rodrigues Parameters \n'
        else:
            out += 'Attitude Type BN    : ERRONEOUS COORDINATE \n'
        out += 'Attitude Coord BR   : '+str(np.around(self.attBR,5))+' \n'
        out += 'Angular Velocity BR : '+str(np.around(self.ohmBR,5))+' rad/s\n'
        out += '\n'
        
        # Print out the current forces enabled
        out += 'FORCES \n'
        out += '------ \n'
        for key in self.forces:
            force_str = str(key)
            while len(force_str) < 20:
                force_str += ' '
            force_str += ': ' + str(self.forces[key]) + ' \n'
            out += force_str
        print(out)
        return out
    
    def _update_states(self):
        '''
        Method that updates the inertial positions, velocities, true anomaly
        (radians), mean argument of latitude (radians), value (float), 
        when ingesting six osculating Keplerian orbit elements.
        Reference: "Satellite Orbits" by Oliver Montenbruck, Chapter 2.2.3.
        '''
        
        GM = self.GM
        a, e, i, w, R, M = self.a, self.e, self.i, self.w, self.R, self.M
        
        # First, let us solve for the eccentric anomaly.
        eccAnom = anomaly.M2E(M,e)
        
        # With the eccentric anomaly, we can solve for position and velocity
        # in the local orbital frame, using the polar equation for an ellipse.
        # Note that the true anomaly would be included in the computation of
        # position and velocity in the perifocal frame below.
        pos_X = a * ( np.cos(eccAnom) - e)
        pos_Y = a * np.sqrt( 1 - e**2 ) * np.sin(eccAnom)
        pos_norm = np.sqrt( pos_X**2 + pos_Y**2 )
        vel_const = np.sqrt( GM * a ) / pos_norm
        vel_X = vel_const * ( -1 * np.sin(eccAnom) )
        vel_Y = vel_const * ( np.sqrt( 1 - e**2 ) * np.cos(eccAnom) )
        
        # The current coordinates are in the local hill frame, and thus 
        # conversion from hill to inertial would be the transpose of HN.
        DCM_HN = rotation.dcmZ(w) @ rotation.dcmX(i) @ rotation.dcmZ(R)
        DCM_NH = np.transpose(DCM_HN)
        
        # With the hill frame, we can now convert it to the ECI frame.
        pos = DCM_NH @ np.array([ pos_X, pos_Y, 0.0 ]).T
        vel = DCM_NH @ np.array([ vel_X, vel_Y, 0.0 ]).T
        nu = np.arctan2( pos_Y, pos_X )
        
        # Update self values.
        self.__dict__['px'] = pos[0]
        self.__dict__['py'] = pos[1]
        self.__dict__['pz'] = pos[2]
        self.__dict__['vx'] = vel[0]
        self.__dict__['vy'] = vel[1]
        self.__dict__['vz'] = vel[2]
        self.__dict__['states'] = [pos[0],pos[1],pos[2],vel[0],vel[1],vel[2]]
        
        # Update the true anomaly, mean anomaly, mean motion and period.
        self.__dict__['n'] = ( GM / (a**3) )**0.5
        self.__dict__['T'] = 2 * PI * ( (a**3) / GM )**0.5
        self.__dict__['nu'] = nu
        self.__dict__['mu'] = ((nu + w) + PI) % (2*PI) - PI
    
    def _update_elements(self):
        '''
        Method that returns six osculating Keplerian orbit elements
        (a, e, i, w, R, M), and ingests in a 1x3 position vector (inertial)
        and a 1x3 velocity vector (inertial). Reference: "Satellite Orbits"
        by Oliver Montenbruck, Chapter 2.2.4.
        '''
        
        GM  = self.GM
        px, py, pz = self.px, self.py, self.pz
        vx, vy, vz = self.vx, self.vy, self.vz
        
        pos = np.array([px, py, pz])
        vel = np.array([vx, vy, vz])
        
        # First, compute the semi-major axis (assuming closed orbit).
        r = np.linalg.norm(pos)
        a = 1 / ( (2/r) - ( ( (np.linalg.norm(vel))**2 ) / GM ) )
        
        # Second, compute the angular momentum vector of the orbit.
        H      = np.cross(pos,vel)
        H_norm = np.linalg.norm(H)
        H_hat  = H / H_norm
        
        # Third, from the normalised angular momentum, derive the inclination.
        i = np.arctan2( math.sqrt( H_hat[0]**2 + H_hat[1]**2 ), H_hat[2] )
        
        # Fourth, from the normalised angular momentum, derive the RAAN.
        R = np.arctan2( H_hat[0], -1*H_hat[1] )
        
        # Fifth, compute the semi-latus rectum.
        p = ( H_norm**2 / GM )
        
        # Sixth, fetch the mean motion.
        n = (GM / (a**3))**0.5
        
        # Seventh, assuming an elliptical orbit, compute the eccentricity.
        e = math.sqrt( 1 - (p/a) )
        
        # Eighth, compute the eccentric anomaly.
        E = np.arctan2( ( (np.dot(pos,vel) / ((a**2)*n)) ), (1 - r/a) )
        
        # Ninth, we can compute the mean anomaly using Kepler's equation.
        M = E - e*math.sin(E)
        
        # Tenth, the argument of latitude is computed.
        U = np.arctan2( pos[2], ( pos[1]*H_hat[0] - pos[0]*H_hat[1] ) )
        
        # Eleventh, the true anomaly is computed.
        nu = np.arctan2( math.sin(E)*math.sqrt(1-e**2), math.cos(E)-e )
        
        # Twelfth, the argument of perigee is computed.
        w = U - nu
        
        # Update the osculating orbital elements.
        self.__dict__['a'] = a
        self.__dict__['e'] = e
        self.__dict__['i'] = i
        self.__dict__['w'] = w
        self.__dict__['R'] = R
        self.__dict__['M'] = M
        self.__dict__['elements'] = [a,e,i,w,R,M]
        
        # Update the true anomaly, mean anomaly, mean motion and period.
        self.__dict__['nu'] = nu
        self.__dict__['mu'] = nu + w
        self.__dict__['n'] = n
        self.__dict__['T'] = 2 * PI / n
    
    # Compute the Hill-Frame transformation matrix.
    def get_hill_frame(self):
        pC = [self.px, self.py, self.pz]
        vC = [self.vx, self.vy, self.vz]
        hC = np.cross(pC, vC) # Angular momentum vector                      
        r_hat = pC / np.linalg.norm(pC) # Local X-axis
        h_hat = hC / np.linalg.norm(hC) # Local Z-axis
        y_hat = np.cross(h_hat, r_hat)  # Local Y-axis
        return np.array([r_hat, y_hat, h_hat])
    
    # Propagate the attitude of the vehicle given a reference and torque.
    # Note that torques and angular velocities must be in the body frame.
    def propagate_attitude(self, dt, torque):
        
        # Get the angular acceleration from Euler's EOM.
        inertia_inverse = np.linalg.inv( self.inertia )
        gyroscopic = np.cross( self.ohmBN, self.inertia @ self.ohmBN )
        wDotBN = inertia_inverse @ ( torque - gyroscopic )
        
        # Check if the coordinate type is a quaternion.
        if self.attBN.strID() == 'QTR' and self.attBR.strID() == 'QTR':
            if self.attBN[0] < 0.0:
                self.attBN.qtr = -1 * self.attBN.qtr # Fix long/short rotation
            qDotBN = self.attBN.get_qtrRate( self.ohmBN )
            self.ohmBN = self.ohmBN + ( dt * wDotBN )
            self.attBN.qtr = self.attBN.qtr + ( dt * qDotBN )
            self.attIntgErr = self.attIntgErr + ( self.attBR.qtr[1:] * dt )
        
        # # Check if the coordinate type is an MRP
        # elif self.attBN.strID() == 'MRP' and self.attBR.strID() == 'MRP':
        #     if self.attBN[0] < 0.0:
        #         self.attBN.qtr = -1 * self.attBN.qtr # Fix long/short rotation
        #     qDotBN = self.attBN.get_qtrRate( self.ohmBN )
        #     self.ohmBN = self.ohmBN + ( dt * wDotBN )
        #     self.attBN.qtr = self.attBN.qtr + ( dt * qDotBN )
        #     self.attIntgErr = self.attIntgErr + ( self.attBR.qtr[1:] * dt )
        else:
            
            # NEED TO IMPLEMENT AN ATTITUDE PROPAGATION FOR CRPs AND MRPs!
            
            print('No quaternion detected in attitude propagation!')
            raise TypeError('MRP and CRP not implemented yet! To be fixed.')
    
    # Two-body Keplerian propagator.
    def propagate_orbit(self, dt):
        self.M = R2D * (self.M + self.n * dt)
        self.epoch += datetime.timedelta( seconds = dt )
        
        
    # # Numerical propagator using RK4.
    # def propagate_perturbed(self, t, step, integrator='RK4'):
    #     time_left = t
    #     if integrator == 'RK4':
    #         while time_left > step:
    #             self.history_ephemeris[self.epoch] = self.states+self.elements
    #             integrate.RK4(self, step)
    #             self.epoch += datetime.timedelta( seconds = step )
    #             time_left -= step
    #         self.history_ephemeris[self.epoch] = self.states+self.elements
    #         integrate.RK4(self, time_left)
    #         self.epoch += datetime.timedelta( seconds = time_left )
    
    def plot_orbit(self):
        return None
    
    def plot_groundtrack(self):
        return None
    