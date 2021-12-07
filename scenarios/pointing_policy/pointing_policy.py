# -*- coding: utf-8 -*-
"""
Reinforcement Learning for Attitude Pointing Path Optimization
in Satellite Formation Flying

Samuel Y. W. Low. Created 28th November 2021, Stanford University.

"""



###############################################################################
###############################################################################
###                                                                         ###
###       Import libraries and change the working directory to root.        ###
###                                                                         ###
###############################################################################
###############################################################################

# Import global libraries
import os
import csv
import math
import copy
import itertools
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

# Close all pre-existing plots.
#plt.close("all")

# Check if we are currently in the root directory of `quadrant`.
while os.getcwd().split("\\")[-1] != "quadrant":
    os.chdir("..") # If not, then move back up one directory.

# Import local source libraries
from source import spacecraft
from source import attitudes
from source import targeting
from source import feedback

action_array = []
time_stamp_array = []

###############################################################################
###############################################################################
###                                                                         ###
###             Import all spacecraft ephemeris and parameters.             ###
###                                                                         ###
###############################################################################
###############################################################################

# Import spacecraft parameters.
sDs = [] # List of deputy spacecraft
deputy_path = os.getcwd() + '\\scenarios\\pointing_policy\\elements.csv'

# Open the CSV file with deputy spacecraft elements computed by QLUSTER.
with open( deputy_path ) as deputy_elements:
    
    deputy_elements_csv = csv.reader( deputy_elements, delimiter=',' )
    
    for row in deputy_elements_csv:
        name = row[0]
        try:
            a    = float(row[1])
            e    = float(row[2])
            i    = float(row[3])
            w    = float(row[4])
            R    = float(row[5])
            M    = float(row[6])
        except:
            pass
        
        # Grab the chief satellite parameters.
        if 'Chief' in name:
            sC = spacecraft.Spacecraft( elements = [a, e, i, w, R, M] )
        
        # Grab the deputy satellites parameters.
        if 'Deputy' in name:
            sD = spacecraft.Spacecraft( elements = [a, e, i, w, R, M] )
            sDs.append( sD )

# `sDs` is a list of spacecraft.Spacecraft() deputy satellites, which
# we can iterate over for solving our Markov Decision Process.



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
aBN = attitudes.QRT( q = [1,0,0,0] )  # Chief body-inertial attitude.
wBN = np.array([0.0, 0.0, 0.0])       # Chief body-inertial angular velocity.

# Initialise arrays for plotting the dynamics.
aBR_array = np.empty((0,4), float)    # Quarternions body-to-reference 
aBN_array = np.empty((0,4), float)    # Quarternions body-to-inertial 
wBR_array = np.empty((0,4), float)    # Angular velocities, body-to-reference
wBN_array = np.empty((0,4), float)    # Angular velocities, body-to-inertial
trq_array = np.empty((0,4), float)    # Control torque applied, inertial frame



###############################################################################
###############################################################################
###                                                                         ###
###  Function that simulates pointing. If used for hallucinating rollouts,  ###
###  then the spacecraft and attitude objects should be deep copies so as   ###
###  to prevent over-writing original mutable objects in the actual sim.    ###
###                                                                         ###
###############################################################################
###############################################################################

def simulate_pointing( dt, ct, sC, sD, Kp, Ki, Kd, aBN, wBN,
                       inertia, torq, intgErr, N ):
    
    print("Performing a rollout for the pointing at Deputy: ", str(N)," \n")
    
    # Initialise a virtual time.
    t = 0.0
    dcmRN = np.eye(3)
    
    # Loop until attitude error converges to zero.
    while True:
    
        # Get the reference target deputy.
        dcmRN, wRN = targeting.reference_deputy( dt, sC, sD, dcmRN )
        
        # Compute the feedback errors and control torques.
        aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QRT( t, dt, ct, 
            Kp, Ki, Kd, aBN, wBN, wRN, dcmRN, inertia, torq, intgErr )
        
        # Propagate both the spacecrafts
        sC.twobody_propagate( dt )
        sD.twobody_propagate( dt )
        
        # Update time
        t += dt
        
        # Successful attitude manoeuvre defined as ten half-lives worth of
        # the average settling time, or when reaching ~ 0.1% in attitude
        # error. Check if the manoeuvre is within this tolerance.
        if sum( np.abs( aBR[1:] ) ) < 0.001 and abs( aBR[0] ) > 0.999:
            break
    
    # Return the final chief and deputy spacecraft and time taken.
    return t, sC, sD

###############################################################################
###############################################################################

def simulate_sun_pointing( dt, ct, sC, Kp, Ki, Kd, aBN, wBN,
                           inertia, torq, intgErr ):
    
    print("Performing a rollout for the sun-pointing manoeuvre. \n")
    
    # Initialise a virtual time.
    t = 0.0
    dcmRN = np.eye(3)
    
    # Loop until attitude error converges to zero.
    while True:
    
        # Get the reference target deputy.
        dcmRN, wRN = targeting.reference_sun()
        
        # Compute the feedback errors and control torques.
        aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QRT( t, dt, ct, 
            Kp, Ki, Kd, aBN, wBN, wRN, dcmRN, inertia, torq, intgErr )
        
        # Propagate both the spacecrafts
        sC.twobody_propagate( dt )
        
        # Update time
        t += dt
        
        # Successful attitude manoeuvre defined as ten half-lives worth of
        # the average settling time, or when reaching ~ 0.1% in attitude
        # error. Check if the manoeuvre is within this tolerance.
        if sum( np.abs( aBR[1:] ) ) < 0.001 and abs( aBR[0] ) > 0.999:
            break
    
    # Return the final chief and deputy spacecraft and time taken.
    return t, sC



###############################################################################
###############################################################################
###                                                                         ###
###     Functions for computing dynamic transitions and rewards for RL      ###
###                                                                         ###
###############################################################################
###############################################################################

# Define the transition matrix function for this scenario.
def transition( n, P, cP ):
    ''' Generates a transition matrix for spacecraft state transitions, that
    depend on the current power levels (or an equivalent conserved constraint)
    in the chief spacecraft.

    Parameters
    ----------
    n : int
        Number of deputy satellites.
    P : float
        Power level of the chief satellite (must be between 0 and 1)
    cP : float
        Power charging urgency constant. Higher values implies a lower
        chance of the spacecraft transiting into charging mode at urgently
        lower power levels.

    Returns
    -------
    T : numpy.ndarray
        Transition matrix of size (N+1) x (N+1), where N is the number of
        deputy satellites and N+1 is the total number of states.
    
    '''
    
    epsilon = ( 1.0 - P ) ** cP
    transition_matrix = ( 1.0 - epsilon ) * np.eye( n+1 )
    transition_matrix[0 ,0] = 1.0
    transition_matrix[1:,0] = epsilon
    return transition_matrix

###############################################################################
###############################################################################

# Define the softmax rewards vector for this scenario.
def reward( P, cP, mu, charge_duration, slew_times ):
    ''' Immediate rewards vector; dynamic rewards that depend on power levels
    and hallucinated slewing times for attitude control manoeuvres. Requires a
    list of attitude manoeuvre slew times, therefore, the agent must first
    perform a single bound hallucination prior to calling the rewards function
    so as to determine the reward distribution.
    
    Parameters
    ----------
    n : int
        Number of deputy satellites.
    P : float
        Power level of the chief satellite (must be between 0 and 1)
    cP : float
        Power charging urgency constant. Higher values implies a lower
        chance of the spacecraft transiting into charging mode at urgently
        lower power levels.
    mu : float
        Hyper-parameter for the soft-max implementation of the rewards.
    charge_duration : float
        Estimated duration required for charging the spacecraft.
    slew_times : list
        A list of length n of all the slew times required for the chief to
        point to a deputy spacecraft.

    Returns
    -------
    rewards_array : numpy.ndarray
        Row vector of length (n+1), where first element corresponds to the 
        reward of charging, and subsequent elements correspond (in order) 
        to the reward of establishing links with the deputies.

    '''
    
    epsilon = ( 1.0 - max(0.0, P) ) ** cP
    total_time = charge_duration + sum( slew_times )
    slew_times = np.array( slew_times )
    Rp = epsilon / math.exp( mu * ( charge_duration / total_time ) )
    reward = [Rp]
    
    for slew_time in slew_times:
        if slew_time == 0.0:
            reward.append(0.0)
        else:
            Rd = (1.0 - epsilon) / math.exp( mu * ( slew_time / total_time ) )
            reward.append(Rd)
    
    return np.array(reward)

###############################################################################
###############################################################################

# Define a function to perform a Bellman update.
def bellman_update( reward, gamma, transition, utility ):
    '''Bellman equation update, for a specific action (bound to the transition
    model based on that action).

    Parameters
    ----------
    reward : numpy.ndarray
        Array of length (n+1), where n = number of deputy satellites.
    gamma : float
        Discount factor.
    transition : numpy.ndarray
        Matrix of length (n+1) x (n+1), where n = number of deputy satellites.
    utility : numpy.ndarray
        Utility, or value array, from the previous step.

    Returns
    -------
    utility_update : numpy.ndarray
        Updated utility, or value array.
    
    '''
    
    utility_update = reward + gamma * ( transition @ utility )
    return utility_update

###############################################################################
###############################################################################
###                                                                         ###
###            Begin dynamics simulation for the random policy.             ###
###                                                                         ###
###############################################################################
###############################################################################

if mode == modes[0]:
    
    # Hard rule-based charging
    charge_cutoff = 0.2
    
    # Select a random spacecraft
    sDi = np.random.randint( 0, len(sDs) )
    sD  = sDs[ sDi ]
    A = sDi + 1
    
    # Initialise a list to record selected spacecraft.
    selected = []
    actions = []
    timeline = []
    
    # Initialise the time counter.
    t = 0.0
    
    while True:
        
        # Get the reference, either sun-pointing or target deputy.
        if A == 0:
            dcmRN, wRN = targeting.reference_sun()
        else:
            dcmRN, wRN = targeting.reference_deputy( dt, sC, sD, dcmRN )
        
        # Compute the feedback errors and control torques.
        aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QRT( t, dt, ct,
                                                                  Kp, Ki, Kd,
                                                                  aBN, wBN,
                                                                  wRN, dcmRN,
                                                                  inertia,
                                                                  torq,
                                                                  intgErr )
        
        # Update the positions and velocities of the spacecrafts
        sC.twobody_propagate( dt )
        for sDz in sDs:
            sDz.twobody_propagate( dt )
        
        # Update the time and power
        t += dt
        power -= pdrain
        
        # Update matrices for plotting.
        aBN_array = np.vstack( [aBN_array, aBN.q] )
        aBR_array = np.vstack( [aBR_array, aBR.q] )
        wBN_array = np.vstack( [wBN_array, np.append(wBN,  norm(wBN))] )
        wBR_array = np.vstack( [wBR_array, np.append(wBR,  norm(wBR))] )
        trq_array = np.vstack( [trq_array, np.append(torq, norm(torq))] )
        
        # Break the while loop if action is sun-pointing and there are no more
        # deputy spacecraft left to link to.
        if A == 0 and len(sDs) == 0:
            print("In sun-pointing mode, but no deputies left. Exiting.")
            print("Time is ", str(t), "\n")
            break
        
        # Check if there is a need to re-select the spacecraft.
        if sum(np.abs(aBR[1:])) < 0.001 and abs(aBR[0]) > 0.999:
            
            # If the spacecraft is in sun-pointing mode, we need to let the
            # simulation continue running until the charging time is complete.
            # Once charging time is complete, then we bring the chief back
            # into hallucination mode where it will begin search again.
            if A == 0:
                
                power = power + pdrain + pcharge
                
                # Pick a deputy once the the battery charges to full.
                if power >= 1.0:
                    
                    actions.append(A)
                    timeline.append(t)
                    
                    while A in actions:
                        sDi = np.random.randint( 0, len(sDs) )
                        sD = sDs [ sDi ]
                        A = sDi + 1
                    print("Established sun pointing mode. Charge finished.")
                    print("Now picking a random deputy ", str(A))
                    print("Time is ", str(t), "\n")
                    
                    
                    torq = np.zeros(3) # Stop all torque for now.
                    dcmRN = np.eye(3) # Reset the reference DCM.
            
            # If the spacecraft is in deputy-pointing mode, we need to drop
            # the current deputy out of our current state space since our
            # mission has been achieved. Then, we need to bring the chief into
            # hallucination mode again so it will be able to pick a new action
            # or a new deputy or possibly even into sun-pointing mode.
            elif A != 0:
                print("Established communications with Deputy ", str(A))
                selected.append( sDi )
                if len(selected) < len(sDs):
                    
                    actions.append(A)
                    timeline.append(t)
                    
                    # Pick a deputy if there is sufficient battery left.
                    if power > charge_cutoff:
                        while A in actions:
                            sDi = np.random.randint( 0, len(sDs) )
                            sD = sDs [ sDi ]
                            A = sDi + 1
                        torq  = np.zeros(3) # Stop all torque for now.
                        dcmRN = np.eye(3) # Reset the reference DCM.
                        print("Sufficient power, picking deputy ", str(A))
                        print("Time is ", str(t), "\n")
                    
                    # Else, go into sun-pointing mode.
                    else:
                        sDi   = None
                        sD    = None
                        A     = 0
                        torq  = np.zeros(3) # Stop all torque for now.
                        dcmRN = np.eye(3) # Reset the reference DCM.
                        print("Insufficient power, going into sun-pointing.")
                        print("Time is ", str(t), "\n")
                else:
                    actions.append(A)
                    timeline.append(t)
                    print("Finished random pointing scenario.")
                    print("Time is ", str(t), "\n")
                    break
            
###############################################################################
###############################################################################

if mode == modes[0]:
    
    plt.figure(1)
    plt.title("Body-to-Inertial Attitude (Quarternion)")
    plt.plot(aBN_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Attitude Coordinate")
    plt.legend(["Q1", "Q2", "Q3", "Q4"])
    plt.grid()
    
    plt.figure(2)
    plt.title("Body-to-Reference Attitude (Quarternion)")
    plt.plot(aBR_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Attitude Coordinate")
    plt.legend(["Q1", "Q2", "Q3", "Q4"])
    plt.grid()
    
    plt.figure(3)
    plt.title("Body-to-Inertial Angular Velocity (rad/s)")
    plt.plot(wBN_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["wx", "wy", "wz", "Norm"])
    plt.grid()
    
    plt.figure(4)
    plt.title("Body-to-Reference Angular Velocity (rad/s)")
    plt.plot(wBR_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["wX", "wY", "wZ", "Norm"])
    plt.grid()
    
    plt.figure(5)
    plt.title("Torque Distribution (N.m)")
    plt.plot(trq_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["TX", "TY", "TZ", "Norm"])
    plt.grid()



###############################################################################
###############################################################################
###                                                                         ###
###              Begin dynamics simulation for greedy search.               ###
###                                                                         ###
###############################################################################
###############################################################################

if mode == modes[1]:
    
    # Flag for hallucinating the dynamics to some depth in greedy search.
    # This flag must be initialised before the while loop (important!).
    hallucination_mode_on = True
    
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
            
            # Initialise the utility vector.
            U = np.zeros( len(sDs) + 1 )
            
            # Make copies
            sC_copy   = copy.deepcopy( sC   )
            sD_copies = copy.deepcopy( sDs  )
            aBN_copy  = copy.deepcopy( aBN  )
            wBN_copy  = copy.deepcopy( wBN  )
            trq_copy  = copy.deepcopy( torq )
                
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
            R = reward( power, cP, cMu, charge_duration, slew_times_array )
            
            # Compute current state transition.
            T = transition( len(sDs), power, cP )
            
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
        aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QRT( t, dt, ct,
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
        aBN_array = np.vstack( [aBN_array, aBN.q] )
        aBR_array = np.vstack( [aBR_array, aBR.q] )
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
                    
###############################################################################
###############################################################################

if mode == modes[1]:
    
    plt.figure(6)
    plt.title("Body-to-Inertial Attitude (Quarternion)")
    plt.plot(aBN_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Attitude Coordinate")
    plt.legend(["Q1", "Q2", "Q3", "Q4"])
    plt.grid()
    
    plt.figure(7)
    plt.title("Body-to-Reference Attitude (Quarternion)")
    plt.plot(aBR_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Attitude Coordinate")
    plt.legend(["Q1", "Q2", "Q3", "Q4"])
    plt.grid()
    
    plt.figure(8)
    plt.title("Body-to-Inertial Angular Velocity (rad/s)")
    plt.plot(wBN_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["wx", "wy", "wz", "Norm"])
    plt.grid()
    
    plt.figure(9)
    plt.title("Body-to-Reference Angular Velocity (rad/s)")
    plt.plot(wBR_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["wX", "wY", "wZ", "Norm"])
    plt.grid()
    
    plt.figure(10)
    plt.title("Torque Distribution (N.m)")
    plt.plot(trq_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["TX", "TY", "TZ", "Norm"])
    plt.grid()



###############################################################################
###############################################################################
###                                                                         ###
###             Begin dynamics simulation for value iteration.              ###
###                                                                         ###
###############################################################################
###############################################################################

if mode == modes[2]:
    
    # Flag for hallucinating the dynamics to some depth in greedy search.
    # This flag must be initialised before the while loop (important!).
    hallucination_mode_on = True
    
    # Initialise the utility vector.
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
            aBN_copy  = copy.deepcopy( aBN  )
            wBN_copy  = copy.deepcopy( wBN  )
            trq_copy  = copy.deepcopy( torq )
                
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
            R = reward( power, cP, cMu, charge_duration, slew_times_array )
            
            # Compute current state transition.
            T = transition( len(sDs), power, cP )
            
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
        aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QRT( t, dt, ct,
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
        aBN_array = np.vstack( [aBN_array, aBN.q] )
        aBR_array = np.vstack( [aBR_array, aBR.q] )
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
                    
###############################################################################
###############################################################################

if mode == modes[2]:
    
    plt.figure(11)
    plt.title("Body-to-Inertial Attitude (Quarternion)")
    plt.plot(aBN_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Attitude Coordinate")
    plt.legend(["Q1", "Q2", "Q3", "Q4"])
    plt.grid()
    
    plt.figure(12)
    plt.title("Body-to-Reference Attitude (Quarternion)")
    plt.plot(aBR_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Attitude Coordinate")
    plt.legend(["Q1", "Q2", "Q3", "Q4"])
    plt.grid()
    
    plt.figure(13)
    plt.title("Body-to-Inertial Angular Velocity (rad/s)")
    plt.plot(wBN_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["wx", "wy", "wz", "Norm"])
    plt.grid()
    
    plt.figure(14)
    plt.title("Body-to-Reference Angular Velocity (rad/s)")
    plt.plot(wBR_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["wX", "wY", "wZ", "Norm"])
    plt.grid()
    
    plt.figure(15)
    plt.title("Torque Distribution (N.m)")
    plt.plot(trq_array)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Omega Component")
    plt.legend(["TX", "TY", "TZ", "Norm"])
    plt.grid()


# ###############################################################################
# ###############################################################################
# ###                                                                         ###
# ###           Begin dynamics simulation for policy optimisation.            ###
# ###                                                                         ###
# ###############################################################################
# ###############################################################################

# # THERE IS A SERIOUS PROBLEM HERE WHERE THE BELLMAN UPDATE IS NO LONGER
# # APPLICABLE IN THIS STUDY. THUS THIS BECOMES BRUTE FORCE FORWARD SEARCH,
# # WITHOUT CONSIDERATION OF THE POWER CHARGING REQUIREMENTS.

# if mode == modes[2]:
    
#     # Flag for hallucinating the dynamics to some depth in forward search.
#     # This flag must be initialised before the while loop (important!).
#     hallucination_mode_on = True
    
#     # Begin a time counter.
#     t = 0.0
    
#     # Begin dynamics and reinforcement learning simulation
#     while True:
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
#         # Begin "hallucinating" in forward search to some depth d. First, we
#         # need to make copies of our existing variables to prevent variables
#         # in the "hallucination" mode from overwriting our actual variables.
#         if hallucination_mode_on:
            
#             print("Entering forward-search mode. \n")
            
#             # Initialise the utility vector.
#             U = np.zeros( len(sDs) + 1 )
            
#             # Generate all possible search trajectories based on the depth.
#             # Action 0 = sun-pointing, anything else refers to a deputy.
#             action_space = list( range( len(sDs) + 1 ) )
#             trajectories = list( itertools.permutations( action_space, depth ))
#             reward_space = []
            
#             # We will have two for loops. First, to enter a trajectory, and
#             # second, to enter each action.
#             for path in trajectories:
                
#                 # A path is a tuple that looks like
#                 # path = (1,5,0,2) for a depth = 4
#                 # total number of paths = N-Permute-D, where N is the size of
#                 # the action space, and D is the depth.
                
#                 # Make copies of the original spacecraft.
#                 sCi   = copy.deepcopy( sC   ) # Single chief
#                 sDsi  = copy.deepcopy( sDs  ) # List of deputies
#                 aBNi  = copy.deepcopy( aBN  ) # Original attitude error
#                 wBNi  = copy.deepcopy( wBN  ) # Original angular velocity
#                 trqi  = copy.deepcopy( torq ) # Original torque
                
#                 # Enter each trajectory.
#                 for action in path:
                
#                 # Compute immediate reward.
#                 R = reward( power, cP, cMu, charge_duration, slew_times_array )
                
#                 # Compute current state transition.
#                 T = transition( len(sDs), power, cP )
                
#                 # Do bellman update
#                 U = bellman_update( R, G, T, U )
                
                
#             # At this point after exiting both for-loops, the Bellman update
#             # should be done. Switch off hallucination mode.
#             hallucination_mode_on = False
            
#             print("Forward search for current state complete! \n")
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
#             # Select the action that maximizes the reward.
#             # A = 0 implies choose the charging state (i.e. sun-pointing mode).
#             # A = i, where i is any number other than 0, implies choosing a 
#             # deputy, which means index (i-1) in the list of deputies `sDs`.
#             A = np.argmax( U )
            
#             # Print out the results of the forward search.
#             if A == 0:
#                 sD = None
#                 print("Entering sun-pointing mode! Power at: " + str(power))
#                 print("Utility yielded from forward search: " )
#                 for a in range( len(U) ):
#                     print(U[a])
#                 print("Expected charging duration: " + str(charge_duration))
#                 print("\n")
            
#             else:
#                 sD = sDs[ A-1 ] # Pick the deputy.
#                 print("Entering deputy pointing mode! Power at: " + str(power))
#                 print("Utility yielded from forward search = " )
#                 for a in range( len(U) ):
#                     print(U[a])
#                 print("Deputy selected = " + str(A))
#                 print("\n")
        
#         # Now, back to non-hallucination mode. This is the actual dynamics.
#         # Get the reference DCM and angular velocity.
#         if A == 0:
#             dcmRN, wRN = targeting.reference_sun()
#         else:
#             dcmRN, wRN = targeting.reference_deputy( dt, sC, sD, dcmRN )
        
#         # Compute the feedback errors and control torques.
#         aBR, wBR, aBN, wBN, torq, intgErr = feedback.control_QRT( t, dt, ct,
#                                                                   Kp, Ki, Kd,
#                                                                   aBN, wBN,
#                                                                   wRN, dcmRN,
#                                                                   inertia,
#                                                                   torq,
#                                                                   intgErr )
        
#         # Update the positions and velocities of the spacecraft
#         sC.twobody_propagate( dt )
#         for sDz in sDs:
#             sDz.twobody_propagate( dt )
        
#         # Update the time and the power level.
#         t += dt
#         power -= pdrain
        
#         # Update matrices for plotting.
#         aBN_array = np.vstack( [aBN_array, aBN.q] )
#         aBR_array = np.vstack( [aBR_array, aBR.q] )
#         wBN_array = np.vstack( [wBN_array, np.append(wBN,  norm(wBN))] )
#         wBR_array = np.vstack( [wBR_array, np.append(wBR,  norm(wBR))] )
#         trq_array = np.vstack( [trq_array, np.append(torq, norm(torq))] )
        
#         # Break the while loop if action is sun-pointing and there are no more
#         # deputy spacecraft left to link to.
#         if A == 0 and len(sDs) == 0:
#             print("In sun-pointing mode, but no deputies left. Exiting.")
#             print("Time is ", str(t), "\n")
#             print("\n")
#             break
        
#         # If the attitude error with the reference has converged, we need to
#         # bring the chief back into hallucination mode.
#         if sum(np.abs(aBR[1:])) < 0.001 and abs(aBR[0]) > 0.999:
            
#             # If the spacecraft is in sun-pointing mode, we need to let the
#             # simulation continue running until the charging time is complete.
#             # Once charging time is complete, then we bring the chief back
#             # into hallucination mode where it will begin forward search again.
#             if A == 0:
                
#                 power = power + pdrain + pcharge
#                 if power >= 1.0:
#                     hallucination_mode_on = True
#                     print("Established sun pointing mode. Charge finished.")
#                     print("Time is ", str(t), "\n")
            
#             # If the spacecraft is in deputy-pointing mode, we need to drop
#             # the current deputy out of our current state space since our
#             # mission has been achieved. Then, we need to bring the chief into
#             # hallucination mode again so it will be able to pick a new action
#             # or a new deputy or possibly even into sun-pointing mode.
#             elif A != 0:
#                 print("Established communications with Deputy ", str(A))
#                 print("Time is ", str(t), "\n")
#                 del sDs[ A-1 ]
#                 if len(sDs) > 0:
#                     hallucination_mode_on = True
#                     torq = np.zeros(3) # Stop all torque for now.
#                     dcmRN = np.eye(3) # Reset the reference DCM.
#                 else:
#                     print("Finished forward-search pointing scenario.")
#                     print("Time is ", str(t), "\n")
#                     break
                    
# ###############################################################################
# ###############################################################################

# if mode == modes[2]:
    
#     plt.figure(11)
#     plt.title("Body-to-Inertial Attitude (Quarternion)")
#     plt.plot(aBN_array)
#     plt.xlabel("Time (seconds)")
#     plt.ylabel("Attitude Coordinate")
#     plt.legend(["Q1", "Q2", "Q3", "Q4"])
#     plt.grid()
    
#     plt.figure(12)
#     plt.title("Body-to-Reference Attitude (Quarternion)")
#     plt.plot(aBR_array)
#     plt.xlabel("Time (seconds)")
#     plt.ylabel("Attitude Coordinate")
#     plt.legend(["Q1", "Q2", "Q3", "Q4"])
#     plt.grid()
    
#     plt.figure(13)
#     plt.title("Body-to-Inertial Angular Velocity (rad/s)")
#     plt.plot(wBN_array)
#     plt.xlabel("Time (seconds)")
#     plt.ylabel("Omega Component")
#     plt.legend(["wx", "wy", "wz", "Norm"])
#     plt.grid()
    
#     plt.figure(14)
#     plt.title("Body-to-Reference Angular Velocity (rad/s)")
#     plt.plot(wBR_array)
#     plt.xlabel("Time (seconds)")
#     plt.ylabel("Omega Component")
#     plt.legend(["wX", "wY", "wZ", "Norm"])
#     plt.grid()
    
#     plt.figure(15)
#     plt.title("Torque Distribution (N.m)")
#     plt.plot(trq_array)
#     plt.xlabel("Time (seconds)")
#     plt.ylabel("Omega Component")
#     plt.legend(["TX", "TY", "TZ", "Norm"])
#     plt.grid()

# ###############################################################################
# ###############################################################################