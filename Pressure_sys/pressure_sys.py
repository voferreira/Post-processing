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

particle = 'Alumina'
scheme = 'CFDEM-Explicit'

if particle == 'Alginate':
    rp = 0.001332 #            Alumina: 0.00154362 Alginate: 0.001332  #m
    rhop = 1028.805  #         Alumina: 3585.892   Alginate: 1028.805 #kg/m^3
    Np = 107960 #               Alumina: 72401      Alginate: 107960    #[-]
    inlet_velocity = 0.0095   #0.0095 #m/s
    eps_RZ = 0.68 #0.68 #0.55
    eps_exp = 0.646 #0.646 #0.55

elif particle == 'Alumina':
    rp = 0.00154362 #            Alumina: 0.00154362 Alginate: 0.001332  #m
    rhop = 3585.892  #         Alumina: 3585.892   Alginate: 1028.805 #kg/m^3
    Np = 72401 #               Alumina: 72401      Alginate: 107960    #[-]
    inlet_velocity = 0.095   #0.0095 #m/s
    eps_RZ = 0.55
    eps_exp = 0.55

else:
    print(f'Particle not identified: {particle}')
    print('Aborting')
    exit()

g = 9.81        #m/s^2
dp = 2*rp       #m
Vnp = Np * 4/3 * pi * rp**3

rhol = 996.778  #kg/m^3
nu = 0.00000084061 #m^2/s
mu = rhol*nu    #Pa
Db = 0.1        #m
Rb = Db/2       #m
Area = pi*Rb**2 #m^2
Hb = 1.1

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
saveFigDir = currentPath.replace('/VTK', '')

try:
    listVTK.remove('walls')
    listVTK.remove('inlet')
    listVTK.remove('outlet')
    listVTK.remove('VTK')
except:
    print('')

#Create a list of time_steps
time_list = [i.replace('.vtk', '').replace('CFD_', '') for i in listVTK]
time_list = [float(i) for i in time_list]
time_list = [round(2*2.5*i*10e-5, 2) for i in time_list]

#Sort list
time_list, listVTK = (list(t) for t in zip(*sorted(zip(time_list, listVTK))))

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
    
    #Create figure for deltaP vs t
    fig0 = plt.figure()
    ax0 = fig0.add_subplot(111)

    #Create figure for deltaP vs t
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    #Create figure to P at 40cm vs t
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    #creating lists to fill with pressure values
    delta_p_simulated_list = []
    P_40cm = []

    for i in range(0, len(listVTK)):
        #Define df_i as current time step df
        exec(f'df_i = df_{i}') 

        #Convert cell data to point data
        df_i = df_i.cell_centers()

        #Convert cell positions and pressure values to pandas
        df_i_cells_position = pd.DataFrame(df_i.points)
        df_i_p     = pd.DataFrame(df_i['p'])

        #Calculate the avarage pressure at the first layer of cells
        p0 = df_i_p[df_i_cells_position[2] == df_i_cells_position.iloc[0, 2]].mean()[0]

        #The same to the last layer
        pf = df_i_p[df_i_cells_position[2] == df_i_cells_position.iloc[-1, 2]].mean()[0]

        #Find the bed total pressure drop:
        delta_p_simulated = (p0-pf)*rhol
        delta_p_simulated_list.append(delta_p_simulated)

        #Selecting data by its position
        position_data = 40 #cm
        df_i_40cm = df_i_p[df_i_cells_position[2] == position_data] #Z position around positionData
        df_i_40cm = df_i_p[df_i_cells_position[2] < position_data] #first layer below the position
        P_40cm.append(df_i_40cm.mean()[0]*rhol)

    #Plot deltaP vs t
    #ax0.plot(time_list, np.repeat(delta_p_Ergun, len(time_list)))
    fig0.suptitle('Alginate: CFDEM-Explicit')
    ax0.plot(time_list, np.repeat(delta_p_Analytical, len(time_list)))
    ax0.plot(time_list, delta_p_simulated_list,'ok')
    ax0.grid()
    ax0.legend(['Analytical', 'Simulated'])
    ax0.set_ylabel(r'$\Delta P \/\ [Pa]$')
    ax0.set_xlabel(r'$time \/\ [sec]$')

    fig1.suptitle('Alginate: CFDEM-Explicit')
    ax1.plot(time_list[2:], np.repeat(delta_p_Analytical, len(time_list[2:])))
    ax1.plot(time_list[2:], delta_p_simulated_list[2:],'ok')
    ax1.grid()
    ax1.legend(['Analytical', 'Simulated'])
    ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
    ax1.set_xlabel(r'$time \/\ [sec]$')

    #plot P at 40 cm as a function of time

    ax2.plot(time_list, P_40cm, '-ok')
    ax2.legend(['Pressure at 40 cm'])
    ax2.grid()
    ax2.set_ylabel(r'$Pressure \/\ at \/\ 40 \/\ cm \/\ [Pa]$')
    ax2.set_xlabel(r'$time \/\ [sec]$')



    fig0.savefig(f'{saveFigDir}/deltaP_t.png')
    fig1.savefig(f'{saveFigDir}/deltaP_t-Zoom.png')
    fig2.savefig(f'{saveFigDir}/P40_t.png')
    print(f'Pressure as func of time plot was saved to {saveFigDir}')

