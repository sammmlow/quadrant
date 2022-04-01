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

# Policy for one-step look ahead.
def lookahead(ai,d):
    hallucination_mode_on = True
    U = np.zeros( len(sDs) + 1 )
    
    # Begin a time counter.
    t = 0.0
    
    # Initialise a list of deputies that have been selected.
    selected = []
    
    # Begin dynamics and reinforcement learning simulation
    while True:
        
        # Begin "hallucinating" in the greedy search to depth = 1. First, we
        # need to make copies of our existing variables to prevent variables
        # in the "hallucination" mode from overwriting our actual variables.
        if hallucination_mode_on:
            
            print("Entering greedy-search mode. \n")
            
            # Make copies
            sC_copy   = copy.deepcopy( sC   )
            sD_copies = copy.deepcopy( sDs  )
            aBN_copy  = copy.deepcopy( sC.aBN  )
            wBN_copy  = copy.deepcopy( sC.wBN  )
            trq_copy  = copy.deepcopy( sC.torque )
                
            # Initialise the slew times for this depth.
            slew_times_array = []
            
            # Hallucinate sun-pointing at a higher time step.
            sCi  = copy.deepcopy( sC_copy  )
            aBNi = copy.deepcopy( aBN_copy )
            wBNi = copy.deepcopy( wBN_copy )
            trqi = copy.deepcopy( trq_copy )
            tC, sCi = simulate_sun_pointing( dt, ct, sCi, Kp, Ki, Kd, 
                                             aBNi, wBNi, inertia, trqi, 
                                             intgErr )
            
            # Estimate the approximate charging duration.
            charge_duration = tC + ( (1.0 - power) / pcharge )
            
            # Hallucinate pointing at a higher time step for all deputies.
            for k in range( len(sD_copies) ):
                
                # Check if the deputy has been selected before.
                if k in selected:
                    slew_times_array.append(0.0)
                    continue
                
                # Else if it has not, then do a hallucination.
                else:
                
                    # The chief has to copied again as it needs to be re-used
                    # whenever a new deputy is simulated.
                    sCi  = copy.deepcopy( sC_copy  )
                    aBNi = copy.deepcopy( aBN_copy )
                    wBNi = copy.deepcopy( wBN_copy )
                    trqi = copy.deepcopy( trq_copy )
                    
                    # Pick a deputy (no need to deep copy)
                    sDi = sD_copies[ k ]
                    
                    # Simulate the pointing, return the slew time, and the
                    # updated chief and deputy.
                    ti, sCi, sDi = simulate_pointing( dt, ct, sCi, sDi,
                                                      Kp, Ki, Kd, aBN, wBN, 
                                                      inertia, torq,
                                                      intgErr, k+1 )
                    slew_times_array.append(ti)
            
            # Compute immediate reward.
            R = reward( power, λ, μ, charge_duration, slew_times_array )
            
            # Compute current state transition.
            T = transition( len(sDs), power, λ )
            
            # Do bellman update
            U = bellman_update( R, G, T, U )
            
            # At this point after exiting both for-loops, the Bellman update
            # should be done. Switch off hallucination mode.
            hallucination_mode_on = False
            
            print("Greedy search for current state complete! \n")
            
            # Select the action that maximizes the reward.
            # A = 0 implies choose the charging state (i.e. sun-pointing mode).
            # A = i, where i is any number other than 0, implies choosing a 
            # deputy, which means index (i-1) in the list of deputies `sDs`.
            A = np.argmax( U )
            
            # Print out the results of the forward search.
            if A == 0:
                sD = None
                print("Entering sun-pointing mode! Power at: " + str(power))
                print("Utility yielded from forward search: " )
                for a in range( len(U) ):
                    print(U[a])
                print("Expected charging duration: " + str(charge_duration))
                print("\n")
            
            else:
                sD = sDs[ A-1 ] # Pick the deputy.
                print("Entering deputy pointing mode! Power at: " + str(power))
                print("Utility yielded from forward search = " )
                for a in range( len(U) ):
                    print(U[a])
                print("Deputy selected = " + str(A))
                print("\n")
        
        # Now, back to non-hallucination mode. This is the actual dynamics.
        # Get the reference DCM and angular velocity.
        if A == 0:
            dcmRN, wRN = targeting.reference_sun()
        else:
            dcmRN, wRN = targeting.reference_deputy( dt, sC, sD, dcmRN )
        
        # Compute the feedback errors and control torques.
        aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QTR( t, dt, ct,
                                                                  Kp, Ki, Kd,
                                                                  aBN, wBN,
                                                                  wRN, dcmRN,
                                                                  inertia,
                                                                  torq,
                                                                  intgErr )
        
        # Update the positions and velocities of the spacecraft
        sC.twobody_propagate( dt )
        for sDz in sDs:
            sDz.twobody_propagate( dt )
        
        # Update the time and the power level.
        t += dt
        power -= pdrain
        
        # Update matrices for plotting.
        aBN_array = np.vstack( [aBN_array, aBN.qtr] )
        aBR_array = np.vstack( [aBR_array, aBR.qtr] )
        wBN_array = np.vstack( [wBN_array, np.append(wBN,  norm(wBN))] )
        wBR_array = np.vstack( [wBR_array, np.append(wBR,  norm(wBR))] )
        trq_array = np.vstack( [trq_array, np.append(torq, norm(torq))] )
        
        # Break the while loop if action is sun-pointing and there are no more
        # deputy spacecraft left to link to.
        if A == 0 and len(sDs) == 0:
            print("In sun-pointing mode, but no deputies left. Exiting.")
            print("Time is ", str(t), "\n")
            print("\n")
            break
        
        # If the attitude error with the reference has converged, we need to
        # bring the chief back into hallucination mode.
        if sum(np.abs(aBR[1:])) < 0.001 and abs(aBR[0]) > 0.999:
            
            # If the spacecraft is in sun-pointing mode, we need to let the
            # simulation continue running until the charging time is complete.
            # Once charging time is complete, then we bring the chief back
            # into hallucination mode where it will begin forward search again.
            if A == 0:
                
                power = power + pdrain + pcharge
                if power >= 1.0:
                    hallucination_mode_on = True
                    print("Established sun pointing mode. Charge finished.")
                    print("Time is ", str(t), "\n")
            
            # If the spacecraft is in deputy-pointing mode, we need to drop
            # the current deputy out of our current state space since our
            # mission has been achieved. Then, we need to bring the chief into
            # hallucination mode again so it will be able to pick a new action
            # or a new deputy or possibly even into sun-pointing mode.
            elif A != 0:
                print("Established communications with Deputy ", str(A))
                print("Time is ", str(t), "\n")
                #del sDs[ A-1 ]
                selected.append( A-1 )
                if len(selected) < len(sDs):
                    hallucination_mode_on = True
                    torq = np.zeros(3) # Stop all torque for now.
                    dcmRN = np.eye(3) # Reset the reference DCM.
                else:
                    print("Finished greedy-search pointing scenario.")
                    print("Time is ", str(t), "\n")
                    break
    



