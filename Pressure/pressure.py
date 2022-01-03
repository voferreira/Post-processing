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

inlet_velocity = 0.0005 #m/s

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

#Defining case Path
currentPath = os.path.dirname(__file__)
print(currentPath)

#Defining time array:
VTKCounter = len(glob.glob1(f'{currentPath}/CFD/VTK/',"*.vtk"))
dt = 500
end_time = VTKCounter*dt
time_array = np.arange(0, end_time, dt)

#Ploting voidfraction as a func of time
for time in time_array:
    #Reading DF from VTK files
    exec(f'df_{time} = pv.read(f\'{currentPath}/CFD/VTK/CFD_{time}.vtk\')')
    
    #Converting voidfraction data to pandas df
    exec(f'df_{time}_voidfraction = pd.DataFrame(df_{time}.cell_data[\'voidfraction\'])')
    
    #printing voidfraction
    #exec(f'print(df_{time}_voidfraction[df_{time}_voidfraction != 1].mean()[0])')
    
    if figure:
        #ploting
        exec(f'plt.scatter({time}*0.0002, df_{time}_voidfraction[df_{time}_voidfraction != 1].mean()[0], c = \'k\')')

if figure:
    plt.show()

#Defining "df_last" as last time step df
exec(f'df_last = df_{time_array[-1]}') 

#Converting cell data to point data
df_last = df_last.cell_centers()

#Converting cell positions and pressure values to pandas
df_last_cells = pd.DataFrame(df_last.points)
df_last_p     = pd.DataFrame(df_last['p'])

#Calculating the avarage pressure at the first layer of cells
p0 = df_last_p[df_last_cells[2] == df_last_cells.iloc[0, 2]].mean()[0]

#The same to the last layer
pf = df_last_p[df_last_cells[2] == df_last_cells.iloc[-1, 2]].mean()[0]

#Finding the bed total pressure drop:
delta_p_simulated = (p0-pf)*rhol

delta_p_Analytical = Analytical(Np, rp, rhol, rhop, g, Area)
delta_p_Ergun  = Ergun(inlet_velocity, mu, rhol, dp)

print('Pressure differences:')
print(f'Simulation = {delta_p_simulated} | Analytical = {delta_p_Analytical} | Ergun = {delta_p_Ergun}')


if plot_P_t:
    for t in time_array:
        #Defining "df_last" as last time step df
        exec(f'df_last = df_{t}') 

        #Converting cell data to point data
        df_last = df_last.cell_centers()

        #Converting cell positions and pressure values to pandas
        df_last_cells = pd.DataFrame(df_last.points)
        df_last_p     = pd.DataFrame(df_last['p'])

        #Calculating the avarage pressure at the first layer of cells
        p0 = df_last_p[df_last_cells[2] == df_last_cells.iloc[0, 2]].mean()[0]

        #The same to the last layer
        pf = df_last_p[df_last_cells[2] == df_last_cells.iloc[-1, 2]].mean()[0]

        #Finding the bed total pressure drop:
        delta_p_simulated = (p0-pf)*rhol

        plt.scatter(t, delta_p_simulated, c='k')
        plt.scatter(t, delta_p_Ergun, c = 'b')
        #plt.scatter(t, delta_p_Analytical, c = 'r')
    plt.legend(['Simulated', 'Ergun', 'Analytical'])
    plt.show()


if figure:
    fig = plt.figure()
    ax = Axes3D(fig)

    #ax.scatter(pd.DataFrame(mesh.points).iloc[:][0], pd.DataFrame(mesh.points).iloc[:][1], pd.DataFrame(mesh.points).iloc[:][2], c = np.array(mesh.point_data['voidfraction']))
    plt.show()

    # Now plot the grid
    #df_last.cell_data["values"] = df_last.cell_data['voidfraction'] 
    #df_last.plot(show_edges=True)