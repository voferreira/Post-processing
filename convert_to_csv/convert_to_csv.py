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
particle = Lethe_pyvista_tools(simulation_path, prm_file_name='liquid_fluidized_bed.prm')

# Read data
particle.read_lethe_to_pyvista('result_particles.pvd', first=-1)


#############################################################################
'''Save particle data'''

# Loop through all results
for i in range(len(particle.list_vtu)):

    # Store results in 'df'
    exec(f'df = particle.df_{i}')

    # Create the pandas dataframe out of the pyvista df
    exec(f'tab = pd.DataFrame(df.points)')

    # Create an array with the column names
    columns = ['X', 'Y', 'Z']

    # Loop through all arrays
    for array in df.array_names:

        # If it is a vector
        if array in ['Velocity', 'Omega', 'FemForce']:
            
            # Loop through all directions
            for j in range(3):
                
                # Append the name of the column
                columns.append(array + '_' + str(j))

                # Generate a pandas dataframe out of the values in the column
                exec(f'aux = pd.DataFrame(df[array][:, {j}])')

                # Concatenate it to the table
                exec(f'tab = pd.concat([tab, aux], axis = 1)')
        else:
            # Append the name of the column
            columns.append(array)

            # Generate a pandas dataframe out of the values in the column
            exec(f'aux = pd.DataFrame(df[array])')

            # Concatenate it to the table
            exec(f'tab = pd.concat([tab, aux], axis = 1)')
    
    # Correct the name of the columns in the pandas dataframe
    tab.columns = columns

    # Export it as a csv
    tab.to_csv(save_path + '/results_' + str(i) + '.csv')


