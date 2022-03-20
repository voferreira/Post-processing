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
'''Simulation properties'''

particle = 'Alumina'
scheme = 'Model type A - Q1-Q1'

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

inlet_velocity = 0.02

#############################################################################

#############################################################################
'''Functions'''

# Analytical solution to pressure drop in fluidized bed
def Analytical(Np, rp, rhol, rhop, g, Area):
    Mp = Np*4/3*pi*rp**3*rhop
    delta_p = Mp*(rhop - rhol)*g/(rhop * Area)
    return delta_p

#Ergun equation to predict pressure drop at fixed beds
def Ergun(inlet_velocity, mu, rhol, Area, dp, Vnp):
    eps_mf = 0.415
    H_mf = Vnp/((1-eps_mf)*Area)
    delta_p = (150*H_mf*inlet_velocity*mu*(1 - eps_mf)**2/(dp**2*eps_mf**3) + 1.75*H_mf*inlet_velocity**2*rhol*(1 - eps_mf)/(dp*eps_mf**3))
    return delta_p
#############################################################################

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath.replace('/output', '')

#Read name of files in .pvd file
files = pd.read_csv(f'{currentPath}/result_.pvd',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
files = files.dropna()
time_list = files['time'].tolist()
list_vtu = files['vtu'].tolist()
list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in list_vtu]

#Estimating analytical and Ergun equation pressure drop
delta_p_Analytical = Analytical(Np, rp, rhol, rhop, g, Area)
delta_p_Ergun  = Ergun(inlet_velocity, mu, rhol, Area, dp, Vnp)

#Define eps_list and voidfraction_list to append value for each time-step
delta_p_simulated_list = []

#Read VTU data
for i in range(0, len(list_vtu)):

    #Read DF from VTU files
    df = pv.read(f'{currentPath}/{list_vtu[i]}')
    
    #Adjust data format
    #df = df.point_data_to_cell_data()

    #Take the first slice of the domain at pos0
    pos0 = [-0.45, 0, 0]
    slice0 = df.slice(normal=[1, 0, 0], origin = pos0)
    p0 = np.mean(slice0['pressure'])

    #Take the first slice of the domain at posf
    posf = [Hb/2-0.02, 0, 0]
    slicef = df.slice(normal=[1, 0, 0], origin = posf)
    pf = np.mean(slicef['pressure'])

    #Store the total pressure drop
    delta_p_simulated_list.append((p0 - pf)*rhol)

#Export the total pressure results as a function of time
csv = pd.DataFrame([time_list, delta_p_simulated_list], index=['time', 'deltaP']).transpose()
csv.to_csv(f'{saveFigDir}/deltaP_t.csv')


#Plot pressure vs time
fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
fig0.suptitle(f'{particle}: {scheme}')
#ax0.plot(time_list, np.repeat(delta_p_Analytical, len(time_list)), label = 'Analytical')
ax0.plot(time_list, np.repeat(delta_p_Ergun, len(time_list)), label = 'Ergun')
ax0.plot(time_list, delta_p_simulated_list,'ok', label = 'Simulated')
ax0.grid()
ax0.legend()
ax0.set_ylabel(r'$\Delta P \/\ [Pa]$')
ax0.set_xlabel(r'$time \/\ [sec]$')
fig0.savefig(f'{saveFigDir}/deltaP_t.png')

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig1.suptitle(f'{particle}: {scheme}')
ax1.plot(time_list[1:], np.repeat(delta_p_Analytical, len(time_list)-1), '--r',label = 'Analytical')
#ax1.plot(time_list[1:], np.repeat(delta_p_Ergun, len(time_list)-1), label = 'Ergun')
ax1.plot(time_list[1:], delta_p_simulated_list[1:],'ok', label = 'Simulation')
ax1.grid()
ax1.legend()
ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
ax1.set_xlabel(r'$time \/\ [sec]$')
fig1.savefig(f'{saveFigDir}/deltaP_t-Zoom.png')

