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
##    Returns reference DCM and OHM of the target, for attitude reference.   ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 23-Nov-2021 23:01 PM (-8 GMT)                            ##
##    Last modified 15-Mar-2022 13:05 AM (-8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import numpy as np
from numpy.linalg import norm as norm
from source import anomaly
from source import rotation



###############################################################################
###############################################################################
###                                                                         ###
###   Pointing DCM & OHM reference as a deputy satellite in the formation.  ###
###                                                                         ###
###############################################################################
###############################################################################



###############################################################################
###############################################################################
###                                                                         ###
###           Pointing DCM & OHM reference as the Sun, defined to           ###
###              be in the inertial frame Y-axis direction.                 ###
###                                                                         ###
###############################################################################
###############################################################################

def reference_sun():
    
    # For the sun-pointing mode, we have defined the sun to always be
    # located when the main spacecraft points in the n2 unit direction
    # vector (and where n3 is the direction vector pointing out of the
    # ecliptic, and n1 x n2 = n3). The dcmSN is defined below.
    
    # Inputs: None. Since the rotation from inertial to sun-pointing is a
    # fixed rotation (not time-varying), external inputs are not needed.
    
    # Output: Reference attitude in 3x3 DCM and reference angular velocity
    # All outputs taken with respect to the inertial frame (i.e., RN frame).
    
    # Sun-pointing from inertial frame is fixed by definition.
    dcmSN = np.array([ [-1,0,0], [0,0,1], [0,1,0] ])
    ohmSN = [ 0.0, 0.0, 0.0 ]
    
    return dcmSN, ohmSN

###############################################################################
###############################################################################
###                                                                         ###
###               Pointing towards Earth in Nadir direction.                ###
###                                                                         ###
###############################################################################
###############################################################################

def reference_hill( sC ):
    p = np.array([ sC.px, sC.py, sC.pz ])
    v = np.array([ sC.vx, sC.vy, sC.vz ])
    dcmHN = sC.get_hill_frame()
    ohmHN = np.cross(p,v) / np.dot(p,p) # N-frame coordinates
    return dcmHN, ohmHN

def reference_nadir( sC ):
    
    # For the nadir-pointing mode, the nadir is always fixed in definition
    # with respect to the hill frame of the main spacecraft. Thus, we need
    # to compute the hill frame first, and then rotate hill to nadir. The 
    # hill frame needs the S/C position and velocity in the inertial frame.
    # Note, the angular velocity OHM_RN = OHM_HN, because the nadir pointing
    # frame is actually a fixed rotation from the standard hill frame.
    
    dcmHN, ohmHN = reference_hill( sC )
    dcmRH = np.array([ [-1,0,0 ], [0,1,0 ], [0,0,-1] ])
    dcmRN = dcmRH @ dcmHN
    ohmRN = ohmHN
    return dcmRN, ohmRN

def reference_deputy( dt, sC, sD ):
    
    # For the inertial to inter-satellite communications frame, denoted as
    # CN, we need to first resolve the three axes of the pointing frame by
    # determining the relative position vectors from the main to the
    # secondary spacecraft. The antenna of the main spacecraft is assumed
    # to be pointing in the positive b1 axis (body-frame) for now.
    
    # Current relative states
    posCi  = np.array([ sC.px, sC.py, sC.pz ])
    posDi  = np.array([ sD.px, sD.py, sD.pz ])
    posCDi = posDi - posCi
    posXi  = np.cross( posCDi, [0,0,1] )
    
    # Back-propagate the ephemeris by one time step.
    sC.propagate_orbit( -1*dt )
    sD.propagate_orbit( -1*dt )
    
    # Previous relative states
    posCf  = np.array([ sC.px, sC.py, sC.pz ])
    posDf  = np.array([ sD.px, sD.py, sD.pz ])
    posCDf = posDf - posCf
    posXf  = np.cross( posCDf, [0,0,1] )
    
    # Compute current inertial to communications pointing frame
    dcmCN_X = posCDi / norm( posCDi )
    dcmCN_Y = posXi / norm( posXi )
    dcmCN_Z = np.cross( dcmCN_X, dcmCN_Y )
    dcmCN = np.array([ dcmCN_X, dcmCN_Y, dcmCN_Z ])
    
    # Compute previous step inertial to communications pointing frame
    dcmCN_Xf = posCDf / norm( posCDf )
    dcmCN_Yf = posXf / norm( posXf )
    dcmCN_Zf = np.cross( dcmCN_Xf, dcmCN_Yf )
    dcmCNf = np.array([ dcmCN_Xf, dcmCN_Yf, dcmCN_Zf ])
    
    # Estimate the rate of change of the DCM numerically.
    dcmCN_dot = (dcmCN - dcmCNf) / dt
    
    # Finally, we need to compute the angular velocity of relative
    # pointing frame, as seen in the inertial frame.
    ohmCN_tilde = -1 * dcmCN.T @ dcmCN_dot
    ohmCN = [ ohmCN_tilde[2,1], ohmCN_tilde[0,2], ohmCN_tilde[1,0] ]
    
    # Forward propagate the ephemeris by one time step.
    sC.propagate_orbit( dt )
    sD.propagate_orbit( dt )
    
    return dcmCN, ohmCN