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
plot_P_t = False
#############################################################################

#############################################################################
'''Importing Libraries'''
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from sklearn import linear_model

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
saveFigDir = currentPath.replace('/VTK', '')

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
time_list = sorted({round(5*i*10e-6, 2) for i in time_list})

#Reading VTK data
for i in range(0, len(listVTK)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{listVTK[i]}\')')

#Define eps_list to append eps to each time-step
eps_list = []

for i in range(0, len(listVTK)):
    #Define "df_last" as last time step df
    exec(f'df = df_{i}') 

    #Convert cell data to point data
    df = df.cell_centers()

    #Convert cell positions and pressure values to pandas
    df_cells = pd.DataFrame(df.points)
    df_p     = pd.DataFrame(df['p'])

    #Calculate the avarage pressure at the first layer of cells
    z0 = df_cells[df_cells[2] == df_cells.iloc[1296, 2]][2].mean()
    p0 = df_p[df_cells[2] == df_cells.iloc[1296, 2]].mean()[0]

    #Create empty lists to fill with z values and pressure as function of z values
    p_z = []
    z_list = []

    #Loop through z values above z0
    for z in df_cells[2][df_cells[2] >= z0].drop_duplicates().tolist():
        z_list.append(z-z0)
        p_z.append((p0-df_p[df_cells[2] == z].mean()[0])*rhol)

    #Apply least-squares to find the porosity
    regr = linear_model.LinearRegression()
    z_list = np.array(z_list).reshape(-1, 1)
    p_z = np.array(p_z).reshape(-1, 1)
    model = regr.fit(z_list, p_z)
    r_sq = model.score(z_list, p_z)

    x = z_list
    y = p_z

    while r_sq < 0.99:
        x = np.delete(x, -1).reshape(-1, 1)
        y = np.delete(y, -1).reshape(-1, 1)
        model = regr.fit(x, y)
        r_sq = model.score(x, y)

    eps = 1-(model.coef_[0][0]/((rhop-rhol)*g))

    eps_list.append(eps)

#Plot pressure vs z
plt.plot(time_list, eps_list, '--ok')
plt.savefig(f'{saveFigDir}/eps_t.png')

exit()

#The same to the last layer
pf = df_p[df_cells[2] == df_cells.iloc[-1, 2]].mean()[0]

#Finde the bed total pressure drop:
delta_p_simulated = (p0-pf)*rhol

#Estimating analytical and Ergun equation pressure drop
delta_p_Analytical = Analytical(Np, rp, rhol, rhop, g, Area)
delta_p_Ergun  = Ergun(inlet_velocity, mu, rhol, dp)

print('Pressure differences:')
print(f'Simulation = {delta_p_simulated}')

if plot_P_t:
    
    #Create figure for deltaP vs t
    fig0 = plt.figure()
    ax0 = fig0.add_subplot(111)

    #Create figure to P at 40cm vs t
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

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
    ax0.plot(time_list, np.repeat(delta_p_Ergun, len(time_list)))
    ax0.plot(time_list, np.repeat(delta_p_Analytical, len(time_list)))
    ax0.plot(time_list, delta_p_simulated_list,'ok')
    ax0.legend(['Ergun', 'Analytical', 'Simulated'])
    ax0.set_ylabel(r'$\Delta P \/\ [Pa]$')
    ax0.set_xlabel(r'$time \/\ [sec]$')

    #plot P at 40 cm as a function of time
    ax1.plot(time_list, P_40cm, '-ok')
    ax1.legend(['Pressure at 40 cm'])
    ax1.set_ylabel(r'$Pressure \/\ at \/\ 40 \/\ cm \/\ [Pa]$')
    ax1.set_xlabel(r'$time \/\ [sec]$')


    saveFigDir = currentPath.replace('/VTK', '')
    fig0.savefig(f'{saveFigDir}/deltaP_t.png')
    fig1.savefig(f'{saveFigDir}/P40_t.png')
    print(f'Pressure as func of time plot was saved to {saveFigDir}')

