#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira
Date: Sep. 10, 2021
"""
#############################################################################

#############################################################################
'''Parameters'''

#############################################################################

#############################################################################
'''Importing Libraries'''
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv

import os
import sys
#############################################################################

#############################################################################
'''Functions'''

#############################################################################

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath.replace('/test', '')

#Define list of VTK files and time list:
list_vtu = os.listdir(currentPath)
list_vtu = [x for x in list_vtu if "pvd" not in x ]
list_vtu = [x for x in list_vtu if "pvtu" not in x ]

#Create a list of time_steps
time_list = [i.replace('.vtu', '').replace('dam-break_VOF.', '').replace('.', '') for i in list_vtu]
time_list = [float(i) for i in time_list]
time_list = [round(i/10000000, 2) for i in time_list]

#Sort list
time_list, list_vtu = (list(t) for t in zip(*sorted(zip(time_list, list_vtu))))

#Set phase_limit to search for maximum x
phase_limit = 0.1

#Create a list to fill with maximum x in which phase > phase_limit
x_list = []

#Read vtu data
for i in range(0, len(list_vtu)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')

    #Select a data to apply the slice   
    exec(f'df = df_{i}')

    #find max 'x' in which phase > 0
    points = pd.DataFrame(df.points[:, 0])
    phase  = pd.DataFrame(df['phase'])

    x_max = max(points[phase[0] > 0.1].values)[0]
    x_list.append(x_max)

print('List of maximum values for x:')
print(x_list)

fig0 = plt.figure()
fig0.suptitle('Maximum x as function of time')
ax0 = fig0.add_subplot(111)
ax0.plot(time_list, x_list, '-ok')
ax0.set_ylabel(r'$X \/\ [m]$')
ax0.set_xlabel(r'$Time \/\ [sec]$')
fig0.savefig(f'{saveFigDir}/xmax_t.png')
plt.show()