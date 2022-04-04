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
##    Optimal Pointing Sequences in Spacecraft Formation Flying              ##
##    using Online Planning with Resource Constraints.                       ##
##                                                                           ##
##    Plotting code for single scenario results. Demonstration scenario      ##
##    results were obtained by running policy scripts for trials = 1.        ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Advised by Professor Mykel J. Kochenderfer                             ##
##    First created 23-Nov-2021. Last modified 04-Apr-2022.                  ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries.
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

# Recorded results for random search in deterministic scenario.
p1_nodes = ['9', '3', 'S', '6', '8', '7', '5', '2',  'S',  '4'  ]
p1_times = [1938,4844,6154,7522,7818,8090,8384,10388,11565,11841]

# Recorded results for greedy (angle) search in deterministic scenario.
p2_nodes = ['5', 'S', '3', '9', 'S', '6', 'S', '7', '2',  'S',  '4',  '8'  ]
p2_times = [2214,3010,3279,5144,5894,7522,8178,8458,10388,11161,11438,11560]

# Recorded results for greedy (angle rate) search in deterministic scenario.
p3_nodes = ['7','9','8','6', '5', '3', 'S', '4','2' ]
p3_times = [294,591,882,1187,1580,2219,2986,4554,4848]

# Recorded results for one step lookahead in deterministic scenario.
p4_nodes = ['7','6','4','2', '8', '9', 'S', 'S', '3', '5' ]
p4_times = [294,585,840,1142,1624,2038,2786,2788,3060,4847]

# Recorded results for forward search (depth 2) in deterministic scenario.
p5_nodes = ['7','6','4','2', '8', 'S', '3', '5', 'S', '9' ]
p5_times = [294,585,840,1142,1624,2289,2565,4847,5689,5961]

# Recorded results for forward search (depth 3) in deterministic scenario.
p6_nodes = ['7','6','4','2', '8', 'S', '3', '5', 'S', '9' ]
p6_times = [294,585,840,1142,1624,2289,2565,4847,5689,5961]

# Put them all in a nested list.
all_nodes = [p1_nodes, p2_nodes, p3_nodes, p4_nodes, p5_nodes, p6_nodes]
all_times = [p1_times, p2_times, p3_times, p4_times, p5_times, p6_times]

# Relative run times for each script for 1 trial.
runtimes = np.array([8.4, 8.4, 8.4, 24, 97, 870])
print('Relative run-times were:', runtimes/min(runtimes))

# Create the figure subplot object.
plt.rcdefaults()
fig, ax = plt.subplots()
ax.xaxis.grid(True,linestyle='dashed')
ax.set_axisbelow(True)
cmap = plt.cm.viridis
policies_names = ('Random', 'Greedy (Angle)', 'Greedy (Rates)',
                  'Lookahead', 'Fwd Search (d=2)', 'Fwd Search (d=3)')
policies_number = np.arange( len(policies_names) )
for i in range(len(all_nodes)):
    nodes, times = all_nodes[i], all_times[i]
    plt.gca().set_prop_cycle( None )
    for j in range(len(nodes)+1):
        c = cmap( j / np.max(len(nodes)+1)) 
        ax.barh( i, width=times[-1*j], height=0.4, color=c )
        ax.scatter( times[-1*j], i, c='white' )
        if nodes[-1*j] != 'S':
            ax.annotate( nodes[-1*j], (times[-1*j], i+0.45))
        else:
            ax.annotate( nodes[-1*j], (times[-1*j], i+0.45), color='red')

ax.set_yticks( policies_number )
ax.set_yticklabels( policies_names )
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Pointing mission execution time (seconds)')
ax.set_title('Pointing execution time (s) for each policy in demo scenario.')
plt.show()

# # Example data
# people = ('Tom', 'Dick', 'Harry', 'Slim', 'Jim')
# y_pos = np.arange(len(people))
# performance = 3 + 10 * np.random.rand(len(people))
# error = np.random.rand(len(people))

# ax.barh(y_pos, performance, xerr=error, align='center')
# ax.set_yticks(y_pos)
# ax.set_yticklabels(people)





