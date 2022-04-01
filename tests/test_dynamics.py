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
##    Test for functional errors and and numerical accuracy for an example   ##
##    Mars mission. Each function represents a specific test scenario.       ##
##    Solutions from Coursera course: Spacecraft Dynamics and Control.       ##
##    https://www.coursera.org/specializations/spacecraft-dynamics-control   ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##                                                                           ##
##    Last updated:  28-Mar-2022 00:01 (GMT)                                 ##
##    First created: 16-Mar-2022 02:31 (GMT)                                 ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries.
import os
import numpy as np

# Move the current directory up until the Quadrant root.
while os.getcwd().split("\\")[-1] != "quadrant":
    os.chdir("..")

# Import source libraries.
from source import spacecraft, attitudes, targeting, feedback

# Initialize constants, time step and PID control gains.
PI = 3.141592653589793
R2D = 180.0/PI # Radians to degrees
D2R = PI/180.0 # Degrees to radians
dt, Kp, Ki, Kd = 1.0, 0.016, 0.0, 0.4



# Function to test the inertial states of spacecraft orbiting Mars.
def test_ephemeris():
    
    # Set the initial conditions.
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sD = spacecraft.Spacecraft( elements=[20424.20,0,0,0,0,250], GM=42828.3)
    
    # Propagate sC by 450s, sD by 1150s.
    sC.propagate_orbit(450)
    sD.propagate_orbit(1150)
    
    # Test for final conditions.
    sC_pos = np.array([ sC.px, sC.py, sC.pz ])
    sC_pos_true = np.array([-669.28508994, 3227.49826592, 1883.18106617])
    assert abs(sum( sC_pos - sC_pos_true )) < 1E-6
    sD_pos = np.array([ sD.px, sD.py, sD.pz ])
    sD_pos_true = np.array([-5399.15037424, -19697.64252078, 0.0])
    assert abs(sum( sD_pos - sD_pos_true )) < 1E-6
    sC_vel = np.array([ sC.vx, sC.vy, sC.vz ])
    sC_vel_true = np.array([-3.2559645, -0.79778654, 0.21011585])
    assert abs(sum( sC_vel - sC_vel_true )) < 1E-6
    sD_vel = np.array([ sD.vx, sD.vy, sD.vz ])
    sD_vel_true = np.array([1.396568, -0.38280117, 0.0])
    assert abs(sum( sD_vel - sD_vel_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test ephemeris calculation in `Spacecraft()` class successful!")



# Function to test the Hill frame DCM of spacecraft orbiting Mars.
def test_frame_hill():
    
    # Set the initial conditions.
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sC.propagate_orbit(300)
    
    # Compute the inertial to reference DCM and angular velocity.
    dcmRN, ohmRN = targeting.reference_hill( sC )
    
    # Test for final conditions.
    dcmRN_true = np.array([[-0.0464774,   0.87414792,  0.48343072],
                           [-0.98417245, -0.12292213,  0.12765086],
                           [ 0.17101007, -0.46984631,  0.8660254 ]])
    ohmRN_true = np.array([ 0.00015131, -0.00041572,  0.00076626])
    assert abs(sum(sum( dcmRN - dcmRN_true ))) < 1E-6
    assert abs(sum( ohmRN - ohmRN_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for Hill frame in `targeting.py` successful!")



# Function to test the Nadir frame DCM of spacecraft orbiting Mars.
def test_frame_nadir():
    
    # Set the initial conditions.
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sC.propagate_orbit(330)
    
    # Compute the inertial to reference DCM and angular velocity.
    dcmRN, ohmRN = targeting.reference_nadir( sC )
    
    # Test for final conditions.
    dcmRN_true = np.array([[ 0.0725817729, -0.8705775347, -0.4866483765],
                           [-0.9825922052, -0.1460794329,  0.1147752484],
                           [-0.1710100717,  0.4698463104, -0.8660254038]])
    ohmRN_true = np.array([ 0.0001513092, -0.0004157185, 0.0007662564])
    assert abs(sum(sum( dcmRN - dcmRN_true ))) < 1E-6
    assert abs(sum( ohmRN - ohmRN_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for nadir frame in `targeting.py` successful!")



# Function to test sun-pointing frame DCM of spacecraft orbiting Mars.
# Note that the assumption is that solar panels 
def test_frame_sun():
    
    # Set the initial conditions.
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sC.propagate_orbit(330)
    
    # Compute the inertial to reference DCM and angular velocity.
    dcmRN, ohmRN = targeting.reference_sun()
    
    # Test for final conditions.
    dcmRN_true = np.array([[-1,0,0],[0,0,1],[0,1,0]])
    ohmRN_true = np.array([0,0,0])
    assert abs(sum(sum( dcmRN - dcmRN_true ))) < 1E-6
    assert abs(sum( ohmRN - ohmRN_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for sun-pointing frame in `targeting.py` successful!")



# Function to test inter-satellite pointing DCM of two spacecraft.
# Within this function, QTR to MRP conversion is also tested.
def test_frame_deputy():
    
    # Set the initial conditions.
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sD = spacecraft.Spacecraft( elements=[20424.20,0,0,0,0,250], GM=42828.3)
    sC.propagate_orbit(330)
    sD.propagate_orbit(330)
    
    # Compute the inertial to reference DCM and angular velocity.
    dcmRN, ohmRN = targeting.reference_deputy( dt, sC, sD )
    
    # Test for final conditions.
    dcmRN_true = np.array([[-0.2654753864,-0.9609281631,-0.0783574157 ],
                           [-0.9638918114, 0.2662941528, 0.0000000000 ],
                           [ 0.0208661216, 0.0755280714,-0.9969253309 ]])
    ohmRN_true = np.array([1.978316704e-05,-5.465493563e-06, 1.913001007e-04])
    assert abs(sum(sum( dcmRN - dcmRN_true ))) < 1E-6
    assert abs(sum( ohmRN - ohmRN_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for comms pointing frame in `targeting.py` successful!")
          


# Function to test the attitude error for Nadir pointing. This function will
# also test quaternion to MRP attitude coordinate conversion.
def test_pointing_initial_all():
    
    # Set the desired control gains to match with Schaub's scenario.
    Kp, Ki, Kd = 0.0055555556, 0.0, 0.1666666667
    
    # Set the initial conditions.
    initial_mrp = attitudes.MRP( mrp = [0.3,-0.4,0.5] )
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sD = spacecraft.Spacecraft( elements=[20424.20,0,0,0,0,250], GM=42828.3)
    sC.attBN = attitudes.QTR( dcm = initial_mrp.dcm )
    sC.ohmBN = np.array([1.00,1.75,-2.20]) * D2R # Radians per second.
    sC.inertia = np.array([[10,0,0],[0,5,0],[0,0,7.5]]) # kg m^2
    
    # Compute and verify initial attitude errors for sun pointing.
    dcmRN, ohmRN = targeting.reference_sun()
    sC = feedback.control( sC, Kp, Ki, Kd, dcmRN, ohmRN )
    attBR = attitudes.MRP( dcm = sC.attBR.dcm )
    ohmBR = sC.ohmBR
    attBR_true = np.array([-0.7754207665,-0.4738682462, 0.0430789315])
    ohmBR_true = np.array([ 0.0174532925, 0.0305432619,-0.0383972435])
    assert abs(sum( attBR - attBR_true )) < 1E-6
    assert abs(sum( ohmBR - ohmBR_true )) < 1E-6
    
    # Compute and verify initial attitude errors for nadir pointing.
    dcmRN, ohmRN = targeting.reference_nadir( sC )
    sC = feedback.control( sC, Kp, Ki, Kd, dcmRN, ohmRN )
    attBR = attitudes.MRP( dcm = sC.attBR.dcm )
    ohmBR = sC.ohmBR
    attBR_true = np.array([ 0.2622652296, 0.5547045658, 0.0394240510])
    ohmBR_true = np.array([ 0.0168488300, 0.0309287900,-0.0389157600])
    assert abs(sum( attBR - attBR_true )) < 1E-6
    assert abs(sum( ohmBR - ohmBR_true )) < 1E-6
    
    # Compute and verify initial attitude errors for comms pointing.
    dcmRN, ohmRN = targeting.reference_deputy( dt, sC, sD )
    sC = feedback.control( sC, Kp, Ki, Kd, dcmRN, ohmRN )
    attBR = attitudes.MRP( dcm = sC.attBR.dcm )
    ohmBR = sC.ohmBR
    attBR_true = np.array([ 0.2123338400, 0.4142451700,-0.0173578900])
    ohmBR_true = np.array([ 0.0172970900, 0.0306574400,-0.0384368700])
    assert abs(sum( attBR - attBR_true )) < 1E-6
    assert abs(sum( ohmBR - ohmBR_true )) < 1E-6

    # Print a statement of success if assertions pass.
    print("Test for initial attitude error computation successful!")



# Function to test the attitude error for sun pointing. This function will
# also test quaternion to MRP attitude coordinate conversion.
def test_pointing_sun():
    
    # Set the control gains.
    Kp, Ki, Kd = 0.0055555556, 0.0, 0.1666666667
    
    # Set the initial conditions.
    initial_mrp = attitudes.MRP( mrp = [0.3,-0.4,0.5] )
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sC.attBN = attitudes.QTR( dcm = initial_mrp.dcm )
    sC.ohmBN = np.array([1.00,1.75,-2.20]) * D2R # Radians per second.
    sC.inertia = np.diag([10.0,5.0,7.5]) # kg m^2
    
    # Propagate both the orbit and attitude for 400 steps.
    for t in range(401):
        dcmRN, ohmRN = targeting.reference_sun()
        sC = feedback.control( sC, Kp, Ki, Kd, dcmRN, ohmRN )
        sC.propagate_orbit( dt )
        sC.propagate_attitude( dt, sC.torque )
        
        # Validate the body to target attitude errors at 4 specific epochs.
        if t == 15:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([  0.26179429, -0.14179886,  0.47278088])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 100:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([  0.13236412,  0.73646432,  0.66667077])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 200:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([ -0.03736946, -0.71970723, -0.54371711])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 400:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([ -0.00592190,  0.68271043,  0.71901198])
            assert abs(sum( mrp - mrp_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for pointing simulation to the sun ran successfully!")


# Function to test the attitude error for nadir pointing. This function will
# also test quaternion to MRP attitude coordinate conversion.
def test_pointing_nadir():
    
    # Set the control gains.
    Kp, Ki, Kd = 0.0055555556, 0.0, 0.1666666667
    
    # Set the initial conditions.
    initial_mrp = attitudes.MRP( mrp = [0.3,-0.4,0.5] )
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sC.attBN = attitudes.QTR( dcm = initial_mrp.dcm )
    sC.ohmBN = np.array([1.00,1.75,-2.20]) * D2R # Radians per second.
    sC.inertia = np.diag([10.0,5.0,7.5]) # kg m^2
    
    # Propagate both the orbit and attitude for 400 steps.
    for t in range(401):
        dcmRN, ohmRN = targeting.reference_nadir( sC )
        sC = feedback.control( sC, Kp, Ki, Kd, dcmRN, ohmRN )
        sC.propagate_orbit( dt )
        sC.propagate_attitude( dt, sC.torque )
        
        # Validate the body to target attitude errors at 4 specific epochs.
        if t == 15:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([ 0.29647480, -0.18401929,  0.44648561])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 100:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([ 0.62873592, -0.23286795,  0.05279609])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 200:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([-0.64347843,  0.54484444,  0.23533678])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 400:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([-0.64767549,  0.55240762,  0.19151409])
            assert abs(sum( mrp - mrp_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for pointing simulation to the nadir ran successfully!")



# Function to test the attitude error for deputy pointing. This function will
# also test quaternion to MRP attitude coordinate conversion.
def test_pointing_deputy():
    
    # Set the control gains.
    Kp, Ki, Kd = 0.0055555556, 0.0, 0.1666666667
    
    # Set the initial conditions.
    initial_mrp = attitudes.MRP( mrp = [0.3,-0.4,0.5] )
    sC = spacecraft.Spacecraft( elements=[3796.19,0,30,0,20,60], GM=42828.3)
    sD = spacecraft.Spacecraft( elements=[20424.20,0,0,0,0,250], GM=42828.3)
    sC.attBN = attitudes.QTR( dcm = initial_mrp.dcm )
    sC.ohmBN = np.array([1.00,1.75,-2.20]) * D2R # Radians per second.
    sC.inertia = np.diag([10.0,5.0,7.5]) # kg m^2
    
    # Propagate both the orbit and attitude for 400 steps.
    for t in range(401):
        dcmRN, ohmRN = targeting.reference_deputy( dt, sC, sD )
        sC = feedback.control( sC, Kp, Ki, Kd, dcmRN, ohmRN )
        sC.propagate_orbit( dt )
        sC.propagate_attitude( dt, sC.torque )
        sD.propagate_orbit( dt )
        
        # Validate the body to target attitude errors at 4 specific epochs.
        if t == 15:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([ 0.29226133, -0.18303740,  0.44918392])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 100:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([ 0.58394631, -0.26934404,  0.16295878])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 200:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([-0.59987114,  0.72322378,  0.03628461])
            assert abs(sum( mrp - mrp_true )) < 1E-6
        if t == 400:
            mrp = attitudes.MRP( dcm = sC.attBN.dcm )
            mrp_true = np.array([-0.58755407,  0.76975744,  0.02030674])
            assert abs(sum( mrp - mrp_true )) < 1E-6
    
    # Print a statement of success if assertions pass.
    print("Test for pointing simulation to the deputy ran successfully!")



# Execute all tests.
if __name__ == '__main__':
    test_ephemeris()
    test_frame_hill()
    test_frame_nadir()
    test_frame_sun()
    test_frame_deputy()
    test_pointing_initial_all()
    test_pointing_sun()
    test_pointing_nadir()
    test_pointing_deputy()
