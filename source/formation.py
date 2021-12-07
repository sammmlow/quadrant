# -*- coding: utf-8 -*-
"""
Created on Sun Nov  7 16:50:41 2021

@author: sammm
"""

import numpy as np
from source import rotation
from source import spacecraft

class Formation():
    
    def __init__(self, sC, sD):
        '''Initialise a formation object, defined by a chief and a deputy.
        
        Parameters
        ----------
        sC : Spacecraft()
            Chief spacecraft, as an instance of the Spacecraft() class from
            the `spacecraft.py` module.
        sD : Spacecraft()
            Deputy spacecraft, as an instance of the Spacecraft() class from
            the `spacecraft.py` module.

        '''
        
        self.sC = sC # Chief
        self.sD = sD # Deputy
        self.update_ROE( sC, sD )
        self.update_RTN( sC, sD )
            
        # except AttributeError:
        #     print("AttributeError: Check if constructor has 02x Spacecraft()?")
        # except TypeError:
        #     print("TypeError: Inputs must be instances of Spacecraft()!")
        # except:
        #     print("Unknown error occurred. Printing constructor args:")
        #     print(sC, sD)

    # Convert Keplerian orbit elements to relative orbit elements exactly.
    def update_ROE(self, sC, sD):
        self.sC = sC # Chief
        self.sD = sD # Deputy
        self.da = (self.sD.a - self.sC.a) / self.sC.a
        self.dL = (self.sD.M - self.sC.M) + (self.sD.w - self.sC.w)
        self.dL = (self.sD.R - self.sC.R) * np.cos(self.sC.i) + self.dL
        self.ix = (self.sD.i - self.sC.i)
        self.iy = (self.sD.R - self.sC.R) * np.sin(self.sC.i)
        self.ex = (self.sD.e * np.cos(self.sD.w))
        self.ex = (self.ex - (self.sC.e * np.cos(self.sC.w)))
        self.ey = (self.sD.e * np.sin(self.sD.w))
        self.ey = (self.ey - (self.sC.e * np.sin(self.sC.w)))
    
    # Convert inertial frame coordinates to RTN frame coordinates exactly.
    def update_RTN(self, sC, sD):
        self.sC = sC # Chief
        self.sD = sD # Deputy
        pC  = np.array([ self.sC.px, self.sC.py, self.sC.pz ])
        pD  = np.array([ self.sD.px, self.sD.py, self.sD.pz ])
        self.compute_hill_dcm() # Updates self.hill_dcm
        self.RTN = self.hill_dcm @ ( pC - pD )
        return self.RTN
        
    # Compute the Hill-Frame 
    def compute_hill_dcm(self):
        pC = [ self.sC.px, self.sC.py, self.sC.pz ] # Position of chief
        vC = [ self.sC.vx, self.sC.vy, self.sC.vz ] # Velocity of chief
        hC = np.cross(pC, vC)                       # Angular momentum vector
        r_hat = pC / np.linalg.norm(pC)             # Local X-axis
        h_hat = hC / np.linalg.norm(hC)             # Local Z-axis
        y_hat = np.cross(h_hat, r_hat)              # Local Y-axis
        self.hill_dcm = np.array([r_hat, y_hat, h_hat])
        return self.hill_dcm
        
        