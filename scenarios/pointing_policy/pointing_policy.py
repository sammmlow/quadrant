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
##    Scenario: Optimal Pointing Sequences in Spacecraft Formation Flying    ##
##    using Online Planning with Resource Constraints                        ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Advised by Professor Mykel J. Kochenderfer                             ##
##    First created 23-Nov-2021 23:01 PM (-8 GMT)                            ##
##    Last modified 15-Mar-2022 13:05 AM (-8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries.
import numpy as np
import matplotlib.pyplot as plt
import os, csv, math, copy, itertools

# Move the current directory up until the Quadrant root.
while os.getcwd().split("\\")[-1] != "quadrant":
    os.chdir("..")

# Import source libraries.
from source import spacecraft, attitudes, targeting, feedback

# Set the policy: 'random', 'greedy', 'lookahead', 'forward', 'mcts'
mode = 'random'

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

# Initialise dynamic and control time step.
dt, ct = 1.0, 1.0

# Initialise reinforcement learning hyper-parameters.
γ = 0.75 # Gamma: Discount factor in Bellman update ( < 1.0)
λ = 2.0  # Lambda: Power charging urgency constant (float)
μ = 5.0  # Mu: Soft max precision parameter for rewards (float)

# Initialise the resource parameters.
P        = 1.0   # Power level, must be between 0 and 1
P_drain  = 0.001 # Power drain rate, between 0 and 1
P_charge = 0.002 # Power charge rate, between 0 and 1

# Initialise attitude control gains
Kp, Ki, Kd = 0.016, 0.0, 0.4 # Initial (Lyapunov) control gains

# Initialise the action space.
actions = ['sun'] + sDs

# Initialise arrays for plotting the dynamics.
aBR_array = np.empty((0,4), float) # Quarternions body-to-reference 
aBN_array = np.empty((0,4), float) # Quarternions body-to-inertial 
wBR_array = np.empty((0,4), float) # Angular velocities, body-to-reference
wBN_array = np.empty((0,4), float) # Angular velocities, body-to-inertial
trq_array = np.empty((0,4), float) # Control torque applied, inertial frame

# Define the transition matrix function for this scenario.
def transition( n, P, λ ):
    # n = number of deputy spacecraft in action space.
    # P = power level of chief spacecraft, between 0 and 1.
    # λ = power charging priority constant (hyper-parameter)
    ε = (1-P)**λ
    transition_matrix = (1-ε) * np.eye(n+1)
    transition_matrix[0 ,0] = 1.0
    transition_matrix[1:,0] = ε
    return transition_matrix

# Define the soft max rewards vector for this scenario.
def reward( P, λ, μ, C, Δ ):
    # P = power level of chief spacecraft, between 0 and 1.
    # λ = power charging priority hyper-parameter.
    # μ = soft-max hyper-parameter to tune reward hardness.
    # C = expected charging duration of chief spacecraft.
    # Δ = array of all expected slew times to each deputy.
    # δ = individual expected slew time to a particular deputy.
    # t = total duration of all actions (charging + slew).
    ε = ( 1 - max(0,P) )**λ
    t = C + sum( Δ )
    Δ = np.array( Δ )
    Rp = ε / math.exp(μ*C/t)
    R = [Rp] # Array of rewards
    for δ in Δ:
        if δ == 0.0:
            R.append(0.0)
        else:
            Rd = (1-ε) / math.exp(μ*δ/t)
            R.append(Rd)
    return np.array(R)

# Define a function to perform a Bellman update.
def bellman( R, γ, T, U1 ):
    # R  = vector of rewards for each action.
    # γ  = discount factor in Bellman equation.
    # T  = transition matrix.
    # U1 = utility from previous iteration.
    U2 = R + γ * ( T @ U1 )
    return U2

# time constant 2I/P where P is angvel gain (D)

# def forward_search(ai,d):
#     if d == 0:
#         return (None,0)
#     Abest, Ubest = None, -10000
    
#     # For all possible actions
#     for a in sDs:
        
#         # Get the immediate reward of the current state-action.
#         reward = reward( P, λ, μ, C, Δ )
        
#         # Recursive depth-first search, returning best action-value pair.
#         Ap, Up = forward_search(a,d-1)
        
#         # Bellman update of utility from child node after exiting recursion.
#         U = R + bellman( R, γ, T, Up )
        
#         # Pick the highest action-value pair
#         if U > Ubest:
#             Abest, Ubest = a, U
    
#     print('Returning best action-value',Abest,Ubest,'at depth',d,'\n')
#     return (Abest, Ubest)



