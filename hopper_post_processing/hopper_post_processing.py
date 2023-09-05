#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira
Date: May 5th, 2022
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
from Lethe_pyvista_tools import *

import sys
#############################################################################

#############################################################################
'''Simulation properties'''

#Take case path as argument
simulation_path = sys.argv[1]

# Name the save path
save_path = simulation_path

# Create the particle object
particle = Lethe_pyvista_tools(simulation_path, prm_file_name='hopper_wiebke.prm')

# Read data
particle.read_lethe_to_pyvista('out.pvd')


#############################################################################

# Position of the outlet
hopper_outlet = 0

# Create a list to store the number of particles below the outlet
number_of_particles_below = []

# Create a list to store the "flow rate" of particles
rate = []

# Loop through all results
for i in range(len(particle.list_vtu)):

    # Store results in 'df'
    exec(f'df = particle.df_{i}')

    vertical_position = df.points[1]

    # Select the data
    vertical_position_below = vertical_position[vertical_position < hopper_outlet]

    # Number of particles above
    number_of_particles_below.append(len(vertical_position_below))

    if i > 0:
        rate.append((number_of_particles_below[i] - number_of_particles_below[i - 1])/particle.prm_dict['time step'])
    
tab = pd.concat([particle.time_list[1:], rate])
tab.columns = ['time', 'rate']
tab.to_csv(save_path + '/results_.csv')


