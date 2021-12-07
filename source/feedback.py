# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 23:25:33 2021

@author: sammm
"""
import numpy as np
from source import attitudes

def control_QRT( t, dt, ct, Kp, Ki, Kd, qBN, wBN, wRN, dcmRN,
                 inertia, prevTorque, intgErr ):

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
    # --> wRN ... ... Angular velocity vector (rad/s) inertial-to-reference
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
    # seen in the B-frame) of wBN and -wRN. Target wRN can be computed in the 
    # "targeting.py" library, before calling this control function. wBN must
    # be updated in each loop, and thus requires the
    # computation of wDotBN for the update. wDotBN can be computed by
    # providing inputs wBN, the satellite inertia I, and the control law
    # U, into Euler's rotational equation of motion. If omegaDot_RN is needed, 
    # the numerical first order gradient will suffice as an estimate.
    #
    # Take note that wRN and omegaDot_RN must be seen in the B-frame!
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
    
    if np.dot( qBN.q, qBN.q ) != 1.0:
        qBN.q = qBN.q / np.linalg.norm( qBN.q )
        
    # Initialise the reference-to-body and reference-to-inertial attitudes.
    dcmBN = qBN.dcm
    dcmBR = dcmBN @ dcmRN.T
    qBR = attitudes.QRT( dcm = dcmBR )
    
    # Second check: Ensure that the quarternions always have norm = 1.
    if np.dot( qBR.q, qBR.q ) != 1.0:
        qBR.q = qBR.q / np.linalg.norm( qBR.q )
        
    # Check if the current quarternion describes the short or long rotation.
    if qBR[0] < 0.0:
        qBR.q = -1 * qBR.q
        
    # The original wRN is seen in the N-frame. Convert to B-frame.
    wRN = dcmBN @ wRN
    
    # Update omega of body-to-reference (B-frame)
    wBR = wBN - wRN
    
    # Trigger the control response at the control time step `ct`.
    if t % ct  == 0.0:
        torque =  ( -1 * Kp * qBR.q[1:] )
        torque -= ( Kd * wBR )
        torque -= ( Ki * Kp * intgErr[1:] )
        
    # Else, return the control torque from the previous loop.
    else:
        torque = prevTorque
        
    # Solve for wDotBN using Euler's rotational equations of motion.
    inertia_inverse = np.linalg.inv( inertia )
    gyroscopic = np.cross( wBN, inertia @ wBN )
    wDotBN = inertia_inverse @ ( torque - gyroscopic )
    
    # Check if the current quarternion describes the short or long rotation.
    if qBN[0] < 0.0:
        qBN.q = -1 * qBN.q
        
    # Use the quarternion rate method inherent in the quarternion class.
    qDotBN = qBN.compute_qrate( wBN )
    
    # Actual numerical integration via first order derivative update.
    wBN = wBN + ( dt * wDotBN )
    qBN.q = qBN.q + ( dt * qDotBN )
    intgErr = intgErr + ( qBR.q * dt )
    
    # Return the attitude and angular velocity errors with respect to the
    # reference target, with respect to the inertial frame, the integration
    # or accumulated attitude error, and the controller torque.
    return qBR, wBR, qBN, wBN, torque, intgErr