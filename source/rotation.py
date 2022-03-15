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
##    Passive Direction Cosine Matrix Tool Box                               ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 20-May-2021 12:50 PM (+8 GMT)                            ##
##    Last modified 20-May-2021 12:50 PM (+8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import math
import numpy as np

###############################################################################
###############################################################################

def dcmX(t):
    '''Generate the direction cosine matrix for an X-axis rotation of angle t.
    
    Parameters
    ----------
    t : float
        Angle theta (t) is the scalar angle (in radians).

    Returns
    -------
    dcm : numpy.ndarray
        Numpy 3x3 direction cosine matrix.
    
    '''
    
    dcm = np.array([[ 1.0,    0.0,         0.0         ],
                    [ 0.0,    math.cos(t), math.sin(t) ],
                    [ 0.0, -1*math.sin(t), math.cos(t) ]])

    return dcm



def dcmY(t):
    '''Generate the direction cosine matrix for an Y-axis rotation of angle t.
    
    Parameters
    ----------
    t : float
        Angle theta (t) is the scalar angle (in radians).

    Returns
    -------
    dcm : numpy.ndarray
        Numpy 3x3 direction cosine matrix.
    
    '''
    
    dcm = np.array([[ math.cos(t), 0.0, -1*math.sin(t) ],
                    [ 0.0,         1.0,    0.0         ],
                    [ math.sin(t), 0.0,    math.cos(t) ]])
    
    return dcm



def dcmZ(t):
    '''Generate the direction cosine matrix for an Z-axis rotation of angle t.
    
    Parameters
    ----------
    t : float
        Angle theta (t) is the scalar angle (in radians).

    Returns
    -------
    dcm : numpy.ndarray
        Numpy 3x3 direction cosine matrix.
    
    '''
    
    dcm = np.array([[    math.cos(t), math.sin(t), 0.0 ],
                    [ -1*math.sin(t), math.cos(t), 0.0 ],
                    [    0.0,         0.0,         1.0 ]])
    
    return dcm