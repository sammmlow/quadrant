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
import os, csv, math
from copy import deepcopy
import matplotlib.pyplot as plt

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
name_index = 0
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
                sC.name = spacecraft.names[ name_index ]
            if 'Deputy' in row[0]: # Grab the deputy SC parameters.
                sD = spacecraft.Spacecraft( elements = [a,e,i,w,R,M] )
                #sD.name = spacecraft.names[ name_index ]
                sD.name = str(name_index+1)
                sDs.append( sD )
            name_index += 1

# Initialise dynamic and control time step.
dt, ct = 1.0, 1.0

# Initialise reinforcement learning hyper-parameters.
γ = 0.85  # Gamma: Discount factor in Bellman update ( < 1.0)
λ = 2.0   # Lambda: Power charging urgency constant (float)
μ = 0.001 # Mu: Soft max precision parameter for rewards (float)

# Initialise the resource parameters.
P        = 1.0    # Power level, must be between 0 and 1
P_drain  = 0.0002 # Power drain rate per second, between 0 and 1
P_charge = 0.001  # Power charge rate per second, between 0 and 1

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


###############################################################################
###############################################################################


# Define the reward for sun pointing.
def reward_sun( P, λ, μ, C ):
    # P = power level of chief spacecraft, between 0 and 1.
    # λ = power charging priority hyper-parameter.
    # μ = soft-max hyper-parameter to tune reward hardness.
    # C = expected duration of sun pointing plus charging.
    α = ( 1 - max(0,P) ) ** λ
    return α / math.exp(μ*C)


# Define the reward for deputy pointing.
def reward_deputy( P, λ, μ, δ ):
    # P = power level of chief spacecraft, between 0 and 1.
    # λ = power charging priority hyper-parameter.
    # μ = soft-max hyper-parameter to tune reward hardness.
    # δ = individual expected slew time to a particular deputy.
    if P <= 0.01:
        return 0.0
    else:
        α = ( 1 - max(0,P) ) ** λ
        return (1-α) / math.exp(μ*δ)


# Define the plotting of rewards.
def plot_rewards():
    plt.close('all')
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
    plt.title('Reward distribution for a sun-pointing state.')
    plt.contourf(p_grid, t_grid, Z_sun_grid, levels=50)
    plt.colorbar()
    plt.xlabel('State: Current power level prior to state (0 to 1)')
    plt.ylabel('Action: Charging and attitude control duration [sec]')
    plt.figure(2)
    plt.title('Reward distribution for a deputy-pointing state.')
    plt.contourf(p_grid, t_grid, Z_sat_grid, levels=50)
    plt.colorbar()
    plt.xlabel('State: Current power level prior to state (0 to 1)')
    plt.ylabel('Action: Attitude control duration [sec]')


# Define the transition probability between states.
def transition( Sf, Si, P, λ ):
    # Sf = final state
    # Si = initial state
    # P  = power level of chief spacecraft, between 0 and 1.
    # λ  = power charging priority constant (hyper-parameter)
    if P <= 0.0:
        α = 1.0
    else:
        α = (1-P)**λ
    if Sf == 'Sun Pointing' and Si == 'Sun Pointing':
        return 0.0
    elif Sf == 'Sun Pointing' and Si != 'Sun Pointing':
        return α
    elif Sf != 'Sun Pointing' and Si != 'Sun Pointing':
        return 1-α
    else:
        return 1.0


# Define the transition matrix used in value iteration.
def transition_matrix( P, λ, n ):
    # P  = power level of chief spacecraft, between 0 and 1.
    # λ  = power charging priority constant (hyper-parameter)
    # n  = number of deputy spacecraft in action space.
    if P <= 0.0:
        α = 1.0
    else:
        α = (1-P)**λ
    transition_matrix = (1-α) * np.eye(n+1)
    transition_matrix[0 ,0] = 1.0
    transition_matrix[1:,0] = α
    return transition_matrix


# Define a function to perform a Bellman update.
def bellman( Sf, Si, P, λ, R, γ, U1 ):
    # R  = vector of rewards for each action.
    # γ  = discount factor in Bellman equation.
    # T  = transition matrix.
    # U1 = utility from previous iteration.
    T = transition( Sf, Si, P, λ )
    U2 = R + γ * T * U1
    return U2


# Define a function hallucinating slew duration to a deputy. This function 
# will also propagate orbits of all deputies in the simulation.
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


# Define a function hallucinating slew and charge duration to the sun. This 
# function will also propagate orbits of all deputies in the simulation.
def point_sun( dt, P, sCi, sDs ):
    δ = 0.0
    while True: # Attitude maneuver to the sun
        dcmRN, ohmRN = targeting.reference_sun()
        sCi = feedback.control( sCi, Kp, Ki, Kd, dcmRN, ohmRN )
        sCi.propagate_orbit( dt )
        sCi.propagate_attitude( dt, sCi.torque )
        for sDii in sDs:
            sDii.propagate_orbit( dt ) # All spacecraft must move!
        δ += dt
        P -= dt * P_drain
        if sum( np.abs( sCi.attBR[1:] ) ) < 0.001:
            if abs( sCi.attBR[0] ) > 0.999:
                break # Break loop once attitude converges to 0.1% error.
    while P < 1.0: # Begin charging process
        dcmRN, ohmRN = targeting.reference_sun()
        sCi = feedback.control( sCi, Kp, Ki, Kd, dcmRN, ohmRN )
        sCi.propagate_orbit( dt )
        sCi.propagate_attitude( dt, sCi.torque )
        for sDii in sDs:
            sDii.propagate_orbit( dt ) # All spacecraft must move!
        δ += dt
        P += dt * P_charge
    return 1.0, δ # Return power level and duration.


###############################################################################
###############################################################################


# Forward search algorithm to solve for best action-value pair.
def forward_search( si, d, Pi ):
    if d == 0:
        return (None,0)
    a_best, U_best = None, -10000
    if si in actions and si != 'Sun Pointing':
        action_index = actions.index(si)
        actions[action_index] = None
    for sf in actions:
        if sf is not None:
            aBN, aBR = deepcopy(sC.attBN), deepcopy(sC.attBR)
            wBN, wBR = deepcopy(sC.ohmBN), deepcopy(sC.ohmBR)
            if sf == 'Sun Pointing':
                Pf, δ = point_sun( 10*dt, Pi, sC, sDs )
                R = reward_sun( Pi, λ, μ, δ )
            else:
                Pf, δ = point_deputy( 10*dt, Pi, sC, sf, sDs )
                R = reward_deputy( Pi, λ, μ, δ )
            Ap, Up = forward_search( sf, d-1, Pf )
            U = bellman( sf, si, Pf, λ, R, γ, Up )
            if U > U_best:
                a_best, U_best = sf, U
            sC.propagate_orbit( -1*δ )
            sC.attBN, sC.attBR = aBN, aBR
            sC.ohmBN, sC.ohmBR = wBN, wBR
            for sDii in sDs:
                sDii.propagate_orbit( -1*δ )
    if si not in actions and si is not None and si != 'Initial':
        actions[action_index] = si
    #print('Depth =',d,'with parent',si,'and best action-value',a_best,U_best)
    return a_best, U_best


###############################################################################
###############################################################################


# Run the main pointing simulation via forward search below.
if __name__ == '__main__':
    Af, Uf, δf, Tf, Df = 'Initial', 0.0, 0.0, 1, 5
    while len(actions) > 1:
        print('Epoch', Tf, 'with current duration', δf, 'and power', P)
        Af, U = forward_search( 'Initial', min( Df, len(actions) ), P )
        print('Completed search at current duration', δf, 'and power', P)
        print('Optimal action-value pair:', Af, U)
        if Af != 'Sun Pointing':
            P, δi = point_deputy( dt, P, sC, Af, sDs )
            actions.remove( Af )
        else:
            P, δi = point_sun( dt, P, sC, sDs ) 
        δf = δf + δi
        Uf = Uf + U
        print('Completed: duration', δf, 'power', P, 'utility', Uf, '\n')
        Tf = Tf + 1
    