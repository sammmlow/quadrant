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
##    Project: Optimal Pointing Sequences in Spacecraft Formation Flying     ##
##             using Online Planning with Resource Constraints               ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Advised by Professor Mykel J. Kochenderfer                             ##
##    First created 23-Nov-2021 23:01 PM (-8 GMT)                            ##
##    Last modified 15-Mar-2022 13:05 AM (-8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

# Import libraries.
import numpy as np
import matplotlib.pyplot as plt
import os, csv, math, copy, itertools
from source import spacecraft, attitudes, targeting, feedback

# Move the current directory up until the Quadrant root.
while os.getcwd().split("\\")[-1] != "quadrant":
    os.chdir("..")
    
# The chief SC is the variable 'sC'. The list of deputies is sDs.
sDs = []

# Import the ephemeris of all spacecraft from the CSV file.
scenario_path = os.getcwd() + '\\scenarios\\pointing_policy'
with open( scenario_path + '\\ephemeris.csv' ) as ephemeris:
    ephemeris_csv = csv.reader( ephemeris, delimiter=',' )
    for row in ephemeris_csv:
        if 'Header' in row:
            continue
        else:
            a,e,i = float(row[1]), float(row[2]), float(row[3])
            w,R,M = float(row[4]), float(row[5]), float(row[6])
            if 'Chief' in row[0]: # Grab the chief SC parameters.
                sC = spacecraft.Spacecraft( elements = [a,e,i,w,R,M] )
            if 'Deputy' in row[0]: # Grab the deputy SC parameters.
                sD = spacecraft.Spacecraft( elements = [a,e,i,w,R,M] )
                sDs.append( sD )

###############################################################################
###############################################################################
###                                                                         ###
###        Initialise all reinforcement-learning related parameters.        ###
###                                                                         ###
###############################################################################
###############################################################################

# State space: [sP, sD1, sD2, ... sDN ], for N number of deputy spacecraft.
# sP refers to the state of charging (sun-pointing mode) and for some i-th
# deputy spacecraft, sDi refers to performing links with that i-th deputy.

# Initialise reinforcement learning parameters for forward search.
G     = 0.75 # Discount factor in Bellman update (float, < 1.0)
cP    = 2.0  # Power charging urgency constant (float)
cMu   = 5.0  # Soft max precision parameter for rewards (float)
depth = 1    # Depth of forward search (integer)



###############################################################################
###############################################################################
###                                                                         ###
###        Initialise all dynamics-related constants and parameters.        ###
###                                                                         ###
###############################################################################
###############################################################################

# Initialise modes.
modes = ['random', 'greedy', 'policy']

# Select the mode.
mode = 'random' # Randomly picks a satellite
mode = 'greedy' # Picks only the immediate reward
mode = 'policy' # Picks the best value, updates utility via Bellman Equation

# Initialise timing parameters.
dt = 1.0                              # Time step in dynamics loop (s)
ct = 1.0                              # Time step in controls loop (s)

# Initialise all constants.
pi  = math.pi                         # 3.141592653589793
GM  = 398600.4418                     # G * Earth Mass (km**3/s**2)
RE  = 6378.140                        # Earth equatorial radius (km)
D2R = math.pi / 180.0                 # Degree to radian convertor

# Initialise the power parameters.
power   = 1.0                         # Power level, must be between 0 and 1
pdrain  = 0.001                       # Power drain rate, between 0 and 1
pcharge = 0.002                       # Power charge rate, between 0 and 1

# Initialise attitude control parameters (quarternions will be used here).
torq       = np.zeros(3)              # Initial control torque.
dcmRN      = np.eye(3)                # Initial reference-to-inertial DCM
intgErr    = np.zeros(4)              # Initial attitude integration error.
inertia    = np.diag([10,10,10])      # Spacecraft principal inertia tensor.
Kp, Ki, Kd = 0.016, 0.0, 0.4          # Initial (Lyapunov) control gains

# Initialise attitudes and rates.
aBN = attitudes.QTR( qtr = [1,0,0,0] )# Chief body-inertial attitude.
wBN = np.array([0.0, 0.0, 0.0])       # Chief body-inertial angular velocity.

# Initialise arrays for plotting the dynamics.
aBR_array = np.empty((0,4), float)    # Quarternions body-to-reference 
aBN_array = np.empty((0,4), float)    # Quarternions body-to-inertial 
wBR_array = np.empty((0,4), float)    # Angular velocities, body-to-reference
wBN_array = np.empty((0,4), float)    # Angular velocities, body-to-inertial
trq_array = np.empty((0,4), float)    # Control torque applied, inertial frame

