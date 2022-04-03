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
import os, csv
import matplotlib.pyplot as plt
from copy import deepcopy

# Move the current directory up until the Quadrant root.
while os.getcwd().split("\\")[-1] != "quadrant":
    os.chdir("..")

# Import source libraries.
from source import spacecraft, attitudes, targeting, feedback, deputy

# Retrieve the ephemeris
sDs = []
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
                sD = spacecraft.Spacecraft(elements=[a,e,i,w,R,M])
                sDs.append( sD )

# Plot the trajectory for one orbit.
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for sDi in sDs:
    sCi = deepcopy( sC )
    posn_array = []
    for t in range(7200):
        hill = sCi.get_hill_frame()
        posn = hill @ (np.array(sCi.states[:3])-np.array(sDi.states[:3]))
        posn_array.append( posn )
        sCi.propagate_orbit(1)
        sDi.propagate_orbit(1)
        if t == 0:
            ax.scatter(posn[0],posn[1],posn[2])
    posn_array = np.array( posn_array )
    ax.plot(posn_array[:,0],posn_array[:,1],posn_array[:,2])
ax.set_xlabel('Cross-Track [km]')
ax.set_ylabel('In-Track [km]')
ax.set_zlabel('Radial [km]')
ax.scatter(0,0,0,'k')