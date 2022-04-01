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
γ = 0.85  # Gamma: Discount factor in Bellman update ( < 1.0)
λ = 2.0   # Lambda: Power charging urgency constant (float)
μ = 0.001 # Mu: Soft max precision parameter for rewards (float)

# Initialise the resource parameters.
P        = 1.0   # Power level, must be between 0 and 1
P_drain  = 0.001 # Power drain rate per second, between 0 and 1
P_charge = 0.002 # Power charge rate per second, between 0 and 1

# Initialise attitude control gains
Kp, Ki, Kd = 0.016, 0.0, 0.4 # Initial (Lyapunov) control gains

# Initialise the action space.
actions = ['Sun Pointing'] + sDs

# Initialise arrays for plotting the dynamics.
aBR_array = np.empty((0,4), float) # Quarternions body-to-reference 
aBN_array = np.empty((0,4), float) # Quarternions body-to-inertial 
wBR_array = np.empty((0,4), float) # Body angular velocities to reference
wBN_array = np.empty((0,4), float) # Body angular velocities to inertial
trq_array = np.empty((0,4), float) # Control torque applied in body frame



# Define the transition probability between states.
def transition( Sf, Si, P, λ ):
    # Sf = final state
    # Si = initial state
    # P  = power level of chief spacecraft, between 0 and 1.
    # λ  = power charging priority constant (hyper-parameter)
    # n  = number of deputy spacecraft in action space.
    α = (1-P)**λ
    # transition_matrix = (1-α) * np.eye(n+1)
    # transition_matrix[0 ,0] = 1.0
    # transition_matrix[1:,0] = α
    if Sf == 'Sun Pointing' and Si == 'Sun Pointing':
        return 1.0
    if Sf == 'Sun Pointing' and Si != 'Sun Pointing':
        return α
    if Sf != 'Sun Pointing' and Si == 'Sun Pointing':
        return 1-α



# Define the reward for sun pointing.
def reward_sun( P, λ, μ, C ):
    # P = power level of chief spacecraft, between 0 and 1.
    # λ = power charging priority hyper-parameter.
    # μ = soft-max hyper-parameter to tune reward hardness.
    # C = expected duration of sun pointing plus charging.
    α = ( 1 - max(0,P) ) ** λ
    return α / math.exp((μ**1.5)*C)



# Define the reward for deputy pointing.
def reward_deputy( P, λ, μ, δ ):
    # P = power level of chief spacecraft, between 0 and 1.
    # λ = power charging priority hyper-parameter.
    # μ = soft-max hyper-parameter to tune reward hardness.
    # δ = individual expected slew time to a particular deputy.
    if P <= 0.01:
        return -1.0
    else:
        α = ( 1 - max(0,P) ) ** λ
        return (1-α) / math.exp(μ*δ)



# Define a function to perform a Bellman update.
def bellman( Sf, Si, P, λ, R, γ, U1 ):
    # R  = vector of rewards for each action.
    # γ  = discount factor in Bellman equation.
    # T  = transition matrix.
    # U1 = utility from previous iteration.
    T = transition( Sf, Si, P, λ )
    U2 = R + γ * T * U1
    return U2



# Define a function hallucinating slew duration to a deputy at 10X time step.
# This function will also update the orbits of all deputies in the simulation.
def point_deputy( dt, P, sCi, sDi, sDs ):
    δ = 0.0
    while True:
        dcmRN, ohmRN = targeting.reference_deputy( dt, sCi, sDi )
        sCi = feedback.control( sCi, Kp, Ki, Kd, dcmRN, ohmRN )
        sCi.propagate_orbit( dt )
        sCi.propagate_attitude( dt, sCi.torque )
        for sDii in sDs:
            sDii.propagate_orbit( dt ) # All spacecraft must move!
        δ += dt
        if sum( np.abs( sCi.attBR[1:] ) ) < 0.001:
            if abs( sCi.attBR[0] ) > 0.999:
                break # Break loop once attitude converges to 0.1% error.
        
    P -= δ * P_drain
    return P, δ # Return power level and duration.



# Define a function hallucinating slew and charge duration to the sun.
def point_sun( dt, P, sCi, sDs ):
    δ = 0.0
    while True:
        dcmRN, ohmRN = targeting.reference_sun()
        sCi = feedback.control( sCi, Kp, Ki, Kd, dcmRN, ohmRN )
        sCi.propagate_orbit( dt )
        sCi.propagate_attitude( dt, sCi.torque )
        for sDii in sDs:
            sDii.propagate_orbit( dt ) # All spacecraft must move!
        δ += dt
        if sum( np.abs( sCi.attBR[1:] ) ) < 0.001:
            if abs( sCi.attBR[0] ) > 0.999:
                break # Break loop once attitude converges to 0.1% error.
    P -= δ * P_drain
    δ += (1-P)/P_charge
    return 1.0, δ # Return power level and duration.



# Forward search algorithm to solve for best action-value pair.
def forward_search( si, d, Pi ):
    # si = initial pointing choice
    # d  = depth of forward search
    # Pi = initial power level
    
    # At the edge of the leaf node, nothing to choose.
    print('\nCurrent action', si, 'at depth =', d)
    if d == 0:
        print('Reached terminal depth, returning to depth',d+1,'\n')
        return (None,0)
    
    # Initialise a best action-value pair.
    a_best, U_best = None, -10000
    
    # Remove any previously visited state in a particular branch.
    if si in actions and si != 'Sun Pointing':
        print('Removing', si, 'for depth =', d-1, 'replaced with None')
        action_index = actions.index(si)
        actions[action_index] = None
        
        # Save the current spacecraft attitude.
        aBN, aBR = sC.attBN, sC.attBR
    
    # Parse through all possible unvisited states in current branch.
    for sf in actions:
        
        # If the current deputy has not been selected in the branch.
        if sf is not None:
            
            # Get the immediate reward of the current state-action.
            # Use a lower-fidelity hallucination at 10x dt for speed.
            if si == 'Sun Pointing':
                Pf, δ = point_sun( 10*dt, Pi, sC ) # Power and time.
                R = reward_sun( Pi, λ, μ, δ )
            if si != 'Sun Pointing':
                Pf, δ = point_deputy( 10*dt, Pi, sC, si ) # Power and time.
                R = reward_deputy( Pi, λ, μ, δ )
            
            # Recursively enter forward search
            Ap, Up = forward_search( sf, d-1, Pf )
            U = ( sf, si, Pf, λ, R, γ, Up )

            if U > U_best:
                a_best, U_best = a, U
     
    # Return the removed action after back propagating up the branch.
    if si not in actions:
        actions[action_index] = si
        
        # At this point, we need to 'travel back in time'.
        sC.propagate_orbit( -1*δ )
        sC.attBN = aBN
        sC.attBR = aBR
        for sDii in sDs:
            sDii.propagate_orbit( -1*δ )
    
    # Return the best action-value pair for this sub-branch.
    return (a_best, U_best)



# Run the main code below.
if __name__ == '__main__':
    
    # Close all plots
    plt.close('all')
    
    # Plot contour plots of the rewards for resource vs mission.
    p_axis = np.arange( 0.02, 1.0, 0.01 )
    t_axis = np.arange( 0.0, 1500.0, 10.0 )
    p_grid, t_grid = np.meshgrid( p_axis, t_axis )
    Z_sun_grid = np.zeros( np.shape(p_grid) )
    Z_sat_grid = np.zeros( np.shape(p_grid) )
    for i in range( len(p_axis) ):
        for j in range( len(t_axis) ):
            p, t = p_axis[i], t_axis[j]
            Z_sun_grid[j][i] = reward_sun( p, λ, μ, t )
            Z_sat_grid[j][i] = reward_deputy( p, λ, μ, t )
    plt.figure(1)
    plt.title('Reward distribution for sun-pointing across resource levels vs slew time.')
    plt.contourf(p_grid, t_grid, Z_sun_grid, levels=50)
    plt.colorbar()
    plt.figure(2)
    plt.title('Reward distribution for deputy-pointing across resource levels vs slew time.')
    plt.contourf(p_grid, t_grid, Z_sat_grid, levels=50)
    plt.colorbar()
    