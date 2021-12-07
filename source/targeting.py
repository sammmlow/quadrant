# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 23:01:28 2021

@author: sammm
"""

import numpy as np
from source import anomaly
from source import rotation

# Returns reference DCM and OHM of the target, for attitude reference.

###############################################################################
###############################################################################
###                                                                         ###
###   Pointing DCM & OHM reference as a deputy satellite in the formation.  ###
###                                                                         ###
###############################################################################
###############################################################################

def reference_deputy( dt, sC, sD, DCM_PREV ):
    
    # For the inertial to inter-satellite communications frame, denoted as
    # CN, we need to first resolve the three axes of the pointing frame by
    # determining the relative position vectors from the main to the
    # secondary spacecraft. The antenna of the main spacecraft is assumed
    # to be pointing in the positive b1 axis (body-frame) for now.
    
    # Inputs:
    # - Position vectors (inertial) of the chief and deputy (p1 and p2)
    # - DCM of deputy pointing reference [CN] from the previous loop
    # - Flag to indicate first epoch or not, for DCM_CN_DOT estimation
    # - Time step used in the dynamics loop
    
    # Output: Reference attitude in 3x3 DCM and reference angular velocity
    # All outputs taken with respect to the inertial frame (i.e., RN frame).
    
    # Relative position vector, chief pointing to the deputy.
    posC = np.array([ sC.px, sC.py, sC.pz ])
    posD = np.array([ sD.px, sD.py, sD.pz ])
    pos_relative = posD - posC
    
    # Cross the relative position with the inertial Z-axis.
    pos_cross = np.cross( pos_relative, [ 0, 0, 1 ] )
    
    # Get the three axes of the comms pointing frame.
    DCM_CN_X = -1*pos_relative / np.linalg.norm( pos_relative )
    DCM_CN_Y = pos_cross / np.linalg.norm( pos_cross )
    DCM_CN_Z = np.cross( DCM_CN_X, DCM_CN_Y )
    
    # Form the actual DCM describing the formation of this frame.
    DCM_CN = np.array([ DCM_CN_X,  DCM_CN_Y,  DCM_CN_Z ])
    
    # Check for the previous DCM. The initialised DCM should either be
    # an identity or zeros-only matrix, so check for that (this should
    # be documented too).
    if np.trace(DCM_PREV) == 3.0 or np.trace(DCM_PREV) == 0.0:
        DCM_PREVIOUS = DCM_CN
    else:
        DCM_PREVIOUS = DCM_PREV
    
    # We need to estimate the rate of change of how the main-to-
    # secondary spacecraft rotation frame, relative to the inertial,
    # changes with time. The only numerical way to do this is to do a
    # linearised interpolation with the DCM_CN of the previous step.
    DCM_CN_DOT = ( DCM_CN - DCM_PREVIOUS ) / dt
    
    # Finally, we need to compute the angular velocity of relative
    # pointing frame, as seen in the inertial frame.
    OHM_CN_TILDE = -1 * DCM_CN.T @ DCM_CN_DOT
    OHM_CN = [ OHM_CN_TILDE[2,1], OHM_CN_TILDE[0,2], OHM_CN_TILDE[1,0] ]
    
    # Now, let us set the correct reference DCM and angular velocity.
    DCM_RN = DCM_CN # R => Reference, N => Inertial
    OHM_RN = OHM_CN # R => Reference, N => Inertial
    
    return DCM_RN, OHM_RN

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
    # ecliptic, and n1 x n2 = n3). The DCM_SN is defined below.
    
    # Inputs: None. Since the rotation from inertial to sun-pointing is a
    # fixed rotation (not time-varying), external inputs are not needed.
    
    # Output: Reference attitude in 3x3 DCM and reference angular velocity
    # All outputs taken with respect to the inertial frame (i.e., RN frame).
    
    # Sun-pointing from inertial frame is fixed by definition.
    DCM_SN = np.array([[ -1,  0,  0 ],
                       [  0,  0,  1 ],
                       [  0,  1,  0 ]])
    
    # Now, let us set the correct reference DCM and angular velocity.
    DCM_RN = DCM_SN;           # R => Reference, N => Inertial
    OHM_RN = [ 0.0, 0.0, 0.0 ] # R => Reference, N => Inertial
    
    return DCM_RN, OHM_RN

###############################################################################
###############################################################################
###                                                                         ###
###               Pointing towards Earth in Nadir direction.                ###
###                                                                         ###
###############################################################################
###############################################################################

def reference_nadir( sC ):
    
    # For the nadir-pointing mode, the nadir is always fixed in definition
    # with respect to the hill frame of the main spacecraft. Thus, we need
    # to compute the hill frame first, and then rotate hill to nadir. The hill
    # frame requires the S/C position and velocity in the inertial frame.
    
    # Inputs: 04x Angular Keplerian elements (deg), position and velocity (m).
    #         - pos1 -> Inertial frame position vector 1x3 (m)
    #         - vel1 -> Inertial frame velocity vector 1x3 (m)
    #         - i    -> Inclination (degrees)
    #         - R    -> Right Angle of Asc Node (degrees)
    #         - w    -> Argument of Perigee (degrees)
    #         - v    -> True Anomaly (degrees)
    # Output: Reference attitude in 3x3 DCM and reference angular velocity
    # All outputs taken with respect to the inertial frame (i.e., RN frame).
    
    # Initialise and construct the hill frame DCM.
    pos = np.array([ sC.px, sC.py, sC.pz ])
    vel = np.array([ sC.vx, sC.vy, sC.vz ])
    DCM_HN = rotation.dcmZ(sC.w) @ rotation.dcmX(sC.i) @ rotation.dcmZ(sC.R)
    
    # Nadir-pointing frame from hill frame is fixed by definition.
    DCM_RH = np.array([[ -1,  0,  0 ],
                       [  0,  1,  0 ],
                       [  0,  0, -1 ]])
    
    # Since the nadir-pointing frame is a fixed rotation away from the
    # hill frame, we will obtain the inertial-to-nadir frame through
    # the matrix products of the DCMs RH and HN below.
    DCM_RN = DCM_RH @ DCM_HN
    
    # Now, we also need to compute the angular velocity of the main
    # satellite, as seen in the N-frame (inertial). By the definition
    # of angular velocity:
    OHM_HN = np.cross( pos, vel ) / np.dot( pos, pos )
    
    # Note, the angular velocity OHM_RN = OHM_HN, because the nadir
    # pointing frame is actually a fixed rotation from the standard
    # hill frame. If you want to prove this by DCM multiplication,
    # then DCM_NH * d(trueAnom)/dt = DCM_NR * d(trueAnom)/dt. The
    # same answer would have been obtained.
    OHM_RN = OHM_HN
    
    # Now, let us set the correct reference DCM and angular velocity.
    # DCM_RN; # R => Reference, N => Inertial
    # OHM_RN; # R => Reference, N => Inertial
    
    return DCM_RN, OHM_RN
