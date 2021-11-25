#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira
Date: Nov. 25, 2021
"""
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
'''Physical properties'''

g = 9.81        #m/s^2
rp = 250*10**(-6)   #m
dp = 2*rp       #m
rhop = 1000 #kg/m^3
Np = 200000     #[-]
Vnp = Np *4/3*pi*rp**3


rhol = 1  #kg/m^3
nu = 0.00001 #m^2/s
mu = rhol*nu    #Pa
Db = 0.003        #m
Rb = Db/2       #m
Area = pi*Rb**2 #m^2
Hb = 0.2

inlet_velocity = 0.2 #m/s
#############################################################################

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath.replace('/output', '')

#Define list of VTU files and time list:
listVTU = os.listdir(currentPath)
listVTU = [x for x in listVTU if "particles" not in x ]
listVTU = [x for x in listVTU if "pvd" not in x ]
listVTU = [x for x in listVTU if "pvtu" not in x ]

#Reorder listVTU
listVTU = sorted(listVTU)

#Create a list of time_steps
time_list = {i.replace('.vtu', '').replace('result_.', '').replace('.', '') for i in listVTU}
time_list = {float(i) for i in time_list}
time_list = sorted({round(i/10000000, 2) for i in time_list})

#Define eps_list and voidfraction_list to append value for each time-step
total_deltaP = []

#Read VTU data
for i in range(0, len(listVTU)):

    #Read DF from VTU files
    df = pv.read(f'{currentPath}/{listVTU[i]}')
    
    #Adjust data format
    df = df.point_data_to_cell_data()

    #Take the first slice of the domain at pos0
    pos0 = [-0.06, 0, 0]
    slice0 = df.slice(normal=[1, 0, 0], origin = pos0)
    p0 = np.mean(slice0['pressure'])

    #Take the first slice of the domain at posf
    posf = [H/2, 0, 0]
    slicef = df.slice(normal=[1, 0, 0], origin = posf)
    pf = np.mean(slicef['pressure'])

    #Store the total pressure drop
    total_deltaP.append((p0 - pf)*rhol)

#Export the total pressure results as a function of time
csv = pd.DataFrame([time_list, total_deltaP], index=['time', 'deltaP']).transpose()
csv.to_csv(f'{saveFigDir}/deltaP_t.csv')


#Plot pressure vs time
fig1, ax1 = plt.subplots(1, 1)
fig1.suptitle(['Title'])
ax1.plot(time_list, total_deltaP, 'ok', ms = 5)
ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
ax1.set_xlabel(r'$Time \/\ [s]$')
fig1.savefig(f'{saveFigDir}/deltaP_t.png')
