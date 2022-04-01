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
##    FILE DESCRIPTION:                                                      ##
##                                                                           ##
##    This file contains the relative orbit propagation function. It takes   ##
##    the 06x chief and deputy orbit Keplerian elements (in km and radians)  ##
##    as input, and outputs the Hill-frame position and velocity vectors as  ##
##    two sets of Nx3 NumPy matrices, where N = total number of samples.     ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 16-May-2021 15:56 AM (+8 GMT)                            ##
##    Last modified 15-Mar-2022 14:03 PM (-8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import numpy as np

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
        
        try:
            self.sC = sC # Chief
            self.sD = sD # Deputy
            self.update_ROE( sC, sD )
            self.update_RTN( sC, sD )
            
        except AttributeError:
            print("AttributeError: Check if constructor has 02x Spacecraft()?")
        except TypeError:
            print("TypeError: Inputs must be instances of Spacecraft()!")
        except:
            print("Unknown error occurred. Printing constructor args:")
            print(sC, sD)

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
        self.get_hill_dcm() # Updates self.hill_dcm
        self.RTN = self.hill_dcm @ ( pD - pC )
    