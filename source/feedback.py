# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 23:25:33 2021

@author: sammm
"""
import numpy as np
from source import attitudes

def control( sC, Kp, Ki, Kd, dcmRN, ohmRN ):

    # ------------------------------------------------------------------------
    # MEMO ON SIMULATING SPACECRAFT ATTITUDES AND EXECUTION OF CONTROL LAW
    # Written by Samuel Low (07-01-2020), DSO National Laboratories.
    # ------------------------------------------------------------------------
    #
    # The inputs of this function are:
    # --> t   ... ... Current time in the dynamics loop
    # --> dt  ... ... Time step in the dynamics loop
    # --> ct  ... ... Time step in the controls software loop
    # --> Kp     ... ... Feedback control proportional gain
    # --> Ki     ... ... Feedback control integral gain
    # --> Kd     ... ... Feedback control derivative gain
    # --> ctrl_mode   ... ... Control law choice (integer)
    # --> I  ... ... Spacecraft 3x3 inertia matrix
    # --> qBN ... ... Attitude error in inertial-to-body frame [BN]
    # --> wBN ... ... Angular velocity vector (rad/s) inertial-to-reference
    # --> ohmRN ... ... Angular velocity vector (rad/s) inertial-to-reference
    # --> dcmRN ... ... Direction cosine matrix inertial-to-reference [RN]
    # --> intgErr   ... Accumulated qBR error integrated over time
    # --> torq  ... Control torque vector from previous loop (Nm)
    #
    # The outputs of this function are the attitude errors and angular rates:
    # --> qBR ... ... Reference-to-Body attitude
    # --> wBR ... ... Reference-to-Body angular velocity (rad/s)
    # --> qBN ... ... Inertial-to-Body attitude
    # --> wBN ... ... Inertial-to-Body angular velocity (rad/s)
    # --> qBN_DOT ... Inertial-to-Body attitude rate
    # --> intgErr   ... Reference-to-Body attitude integration error
    # --> torque  ... Control torque response (Nm) in current feedback
    #
    # ------------------------------------------------------------------------
    #
    # To understand the feedback control logic in this function, which is to
    # be embedded in the numerical integrator loop, let us work backwards.
    #
    # In general for tracking, the reference-to-body outputs qBR and wBR
    # are to be updated in the feedback. This means qDotBR is required, which
    # can be computed using qBR and wBR using the differential kinematic
    # equation at that instant. wBR can be obtained via vector addition (as
    # seen in the B-frame) of wBN and -ohmRN. Target ohmRN can be computed in the 
    # "targeting.py" library, before calling this control function. wBN must
    # be updated in each loop, and thus requires the
    # computation of wDotBN for the update. wDotBN can be computed by
    # providing inputs wBN, the satellite inertia I, and the control law
    # U, into Euler's rotational equation of motion. If omegaDot_RN is needed, 
    # the numerical first order gradient will suffice as an estimate.
    #
    # Take note that ohmRN and omegaDot_RN must be seen in the B-frame!
    #
    # Keep in mind also that if we are working with MRPs, we should switch to
    # the appropriate shadow and non-shadow set values by checking if the MRP
    # norm exceeds 1. This is not necessary with CRPs, since CRP attitude
    # coordinates are unique (i.e., the shadow set = non-shadow set). In the
    # case of quarternions, it is necessary to check that the Beta-0 term is
    # greater than or equal to 0, in order to ensure that the short rotation
    # is called, instead of the long rotation (think: the -359 degree case).
    #
    # Now, we understand what are the ingredients we need to solve for qBR
    # in each loop, we can update the respective derivatives, retrieve the
    # appropriate reference attitudes, repeat again in the next loop.
    #
    # --------------------------------------------------------------------
    
    # First check: Ensure that the quarternions always have norm = 1.
    # Although quarternions are represented as NumPy arrays, you cannot
    # perform vector addition as quarternions are not like vectors where the
    # law of commutativity and the law of parallelograms is adhered to. Thus,
    # in order to perform vector operations on them, you must call the
    # attribute `q` to draw out the quarternion as a NumPy array.
    
    # If spacecraft attitude coordinates are in quaternions...
    if sC.attBN.strID() == 'QTR' and sC.attBR.strID() == 'QTR':
        
        # Check if BN quaternion describes the short or long rotation.
        if sC.attBN[0] < 0.0:
            sC.attBN.qtr = -1 * sC.attBN.qtr
            
        # Re-normalise BN quaternions if norm is not 1.
        if np.dot( sC.attBN.qtr, sC.attBN.qtr ) != 1.0:
            sC.attBN.normalise()
        
        # Initialise reference-to-body and reference-to-inertial attitudes.
        dcmBN = sC.attBN.dcm
        dcmBR = dcmBN @ dcmRN.T
        sC.attBR = attitudes.QTR( dcm = dcmBR )
        
        # Check if BR quaternion describes the short or long rotation.
        if sC.attBR[0] < 0.0:
            sC.attBR.qtr = -1 * sC.attBR.qtr
            
        # Re-normalise BR quaternions if norm is not 1.
        if np.dot( sC.attBR.qtr, sC.attBR.qtr ) != 1.0:
            sC.attBR.normalise()
        
        # The original ohmRN is seen in the N-frame. Convert to B-frame.
        ohmRN = dcmBN @ ohmRN
        
        # Update angular velocity of body-to-reference (in body frame)
        sC.ohmBR = sC.ohmBN - ohmRN
        
        # Compute the feedback control torque and apply it on spacecraft.
        # Note that this only changes the attribute of the torque on the
        # spacecraft, to meaningfully propagate it, the user should run
        # the "propagate_attitude_stepwise( dt, torque )" method of the SC.
        torque =  ( -1 * Kp * sC.attBR.qtr[1:] )
        torque -= ( Kd * sC.ohmBR )
        torque -= ( Ki * Kp * sC.attIntgErr )
        sC.torque = torque
        return sC
    
    else:
        return False
        
    