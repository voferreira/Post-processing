#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira
Date: Sep. 10, 2021
"""
#############################################################################

#############################################################################
'''Parameters'''

figure = False #True if plot
plot_P_t = True
#############################################################################

#############################################################################
'''Importing Libraries'''
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv

import glob

import os
import sys
#############################################################################

#############################################################################
'''Simulation properties'''

g = 9.81        #m/s^2
rp = 0.001332   #m
dp = 2*rp       #m
rhop = 1028.805 #kg/m^3
Np = 107960     #[-]


rhol = 996.778  #kg/m^3
nu = 0.00000084061 #m^2/s
mu = rhol*nu    #Pa
Db = 0.1        #m
Rb = Db/2       #m
Area = pi*Rb**2 #m^2

inlet_velocity = 0.01 #m/s

#############################################################################

#############################################################################
'''Functions'''

# Analytical solution to pressure drop in fluidized bed
def Analytical(Np, rp, rhol, rhop, g, Area):
    Mp = Np*4/3*pi*rp**3*rhop
    delta_p = Mp*(rhop - rhol)*g/(rhop * Area)
    return delta_p

#Ergun equation to predict pressure drop at fixed beds
def Ergun(inlet_velocity, mu, rhol, dp):
    eps_mf = 0.415
    delta_p = 0.23875*(150*inlet_velocity*mu*(1 - eps_mf)**2/(dp**2*eps_mf**3) + 1.75*inlet_velocity**2*rhol*(1 - eps_mf)/(dp*eps_mf**3))
    return delta_p
#############################################################################

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]

#Define list of VTK files and time list:
listVTK = os.listdir(currentPath)
listVTK.remove('walls')
listVTK.remove('inlet')
listVTK.remove('outlet')
try:
    listVTK.remove('VTK')
except:
    print('')

#Reorder listVTK
listVTK = sorted(listVTK)

#Create a list of time_steps
time_list = {i.replace('.vtk', '').replace('CFD_', '') for i in listVTK}
time_list = {float(i) for i in time_list}
time_list = sorted({5*i*10e-4 for i in time_list})

#Reading VTK data
for i in range(0, len(listVTK)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{listVTK[i]}\')')

#Define "df_last" as last time step df
exec(f'df_last = df_{len(listVTK) - 1}') 

#Convert cell data to point data
df_last = df_last.cell_centers()

#Convert cell positions and pressure values to pandas
df_last_cells = pd.DataFrame(df_last.points)
df_last_p     = pd.DataFrame(df_last['p'])

#Calculate the avarage pressure at the first layer of cells
p0 = df_last_p[df_last_cells[2] == df_last_cells.iloc[1296, 2]].mean()[0]

#The same to the last layer
pf = df_last_p[df_last_cells[2] == df_last_cells.iloc[-1, 2]].mean()[0]

#Finde the bed total pressure drop:
delta_p_simulated = (p0-pf)*rhol

#Estimating analytical and Ergun equation pressure drop
delta_p_Analytical = Analytical(Np, rp, rhol, rhop, g, Area)
delta_p_Ergun  = Ergun(inlet_velocity, mu, rhol, dp)

print('Pressure differences:')
print(f'Simulation = {delta_p_simulated} | Analytical = {delta_p_Analytical} | Ergun = {delta_p_Ergun}')

if plot_P_t:
    for i in range(0, len(listVTK)):
        #Define "df_last" as last time step df
        exec(f'df_last = df_{i}') 

        #Convert cell data to point data
        df_last = df_last.cell_centers()

        #Convert cell positions and pressure values to pandas
        df_last_cells = pd.DataFrame(df_last.points)
        df_last_p     = pd.DataFrame(df_last['p'])

        #Calculate the avarage pressure at the first layer of cells
        p0 = df_last_p[df_last_cells[2] == df_last_cells.iloc[0, 2]].mean()[0]

        #The same to the last layer
        pf = df_last_p[df_last_cells[2] == df_last_cells.iloc[-1, 2]].mean()[0]

        #Find the bed total pressure drop:
        delta_p_simulated = (p0-pf)*rhol

        #Plot results against estimations
        plt.scatter(time_list[i]/100, delta_p_Ergun, c = 'b')
        plt.scatter(time_list[i]/100, delta_p_Analytical, c = 'r')
        plt.scatter(time_list[i]/100, delta_p_simulated, c='k')
        
    plt.legend(['Ergun', 'Analytical', 'Simulated'])
    saveFigDir = currentPath.replace('/VTK', '')
    plt.savefig(f'{saveFigDir}/deltaP_t.png')
    print(f'Pressure as func of time plot was saved to {currentPath}')
exit()

"""
    #Selecting data by its position
    positionData = 40#cm
    exec(f'df_{i}_40cm = df_{i}[round(100*df_{i}[2]) == positionData]') #Z position around positionData
    exec(f'df_{i}_40cm = df_{i}[100*df_{i}[2] < positionData]') #first layer below the position

    plt.scatter(time_list[i]/100, delta_p_simulated, c='k')
    plt.savefig(f'{saveFigDir}/deltaP_t.png')



if figure:
    fig = plt.figure()
    ax = Axes3D(fig)

    #ax.scatter(pd.DataFrame(mesh.points).iloc[:][0], pd.DataFrame(mesh.points).iloc[:][1], pd.DataFrame(mesh.points).iloc[:][2], c = np.array(mesh.point_data['voidfraction']))
    plt.show()

    # Now plot the grid
    #df_last.cell_data["values"] = df_last.cell_data['voidfraction'] 
    #df_last.plot(show_edges=True)"""