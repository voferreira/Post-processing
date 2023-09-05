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
write_to_excel = True
#############################################################################

#############################################################################
'''Importing Libraries'''
from math import e, pi
from re import U
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from tqdm import tqdm
from Lethe_pyvista_tools import *
from openpyxl import load_workbook

import os
import sys
#############################################################################

#############################################################################
'''Simulation properties'''

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath + '/Velocity'

particle_name = sys.argv[2]

fluid = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')
fluid.read_lethe_to_pyvista('result_.pvd', first = 1)
#fluid.read_lethe_to_pyvista_parallel('result_.pvd', last = 200, processors=2)

particle = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')

model_type = fluid.prm_dict['vans model'].replace('model', '')
velocity_order = str(int(fluid.prm_dict['velocity order']))
pressure_order = str(int(fluid.prm_dict['pressure order']))

boundary_condition = fluid.prm_dict['bc 0']

scheme = f'Model type {model_type} - Q{velocity_order}-Q{pressure_order}'

dp = particle.prm_dict['diameter']
rp = dp/2
rhop = particle.prm_dict['density particles']
Np = particle.prm_dict['number']
inlet_velocity = fluid.prm_dict['u']

g = abs(particle.prm_dict['gx'])       #m/s^2
Vnp = Np * 4/3 * pi * rp**3

rhol = fluid.prm_dict['density']  #kg/m^3
nu = fluid.prm_dict['kinematic viscosity'] #m^2/s
mu = rhol*nu    #Pa
Db = 0.1        #m
Rb = Db/2       #m
Area = pi*Rb**2 #m^2
Hb = 1.1

#############################################################################

#############################################################################

#Create a folder to recieve velocity results:
if os.path.isdir(currentPath + '/Velocity') == False:
    os.mkdir(currentPath + '/Velocity')

pbar = tqdm(total = len(fluid.time_list)-1, desc="Processing data")

#Create an empty list to store the normalized velocity
normalized_velocity_list = []
normalized_voidfraction_list = []

for i in range(len(fluid.time_list)):
    #Define "df" as last time step df
    exec(f'df = fluid.df_{i}') 
    
    #Take the first slice of the domain at pos
    position = [0.05, 0, 0]

    #Slice df in the middle of the bed
    slice = df.slice(position)

    #Create a radius vector
    slice['R'] = (slice.points[:, 1]**2 + slice.points[:, 2]**2)**(0.5)

    #Create a set with unrepeated radius
    radius = set(slice['R'])

    #Create an upwards velocity vector:
    slice['velocity_X'] = slice['velocity'][:, 0]

    #Calculate and store the normalized velocity
    normalized_velocity = np.mean(slice['velocity_X'][slice['R'] == min(slice['R'])])/np.mean(slice['velocity_X'][slice['R'] == max(slice['R'])])
    normalized_velocity_list.append(normalized_velocity)

    normalized_voidfraction = np.mean(slice['void_fraction'][slice['R'] == min(slice['R'])])/np.mean(slice['void_fraction'][slice['R'] == max(slice['R'])])
    normalized_voidfraction_list.append(normalized_voidfraction)

    #Create an empty list to store the velocity as a radial function
    velocity_radius_list = []
    voidfraction_radius_list = []

    #Create a set of radius
    radius_set = list(set(slice['R']))
    radius_set = np.sort(radius_set)

    for r in radius_set:
        velocity_radius_list.append(np.mean(slice['velocity_X'][slice['R'] == r]))

    '''#Store velocity as a radial function
    fig2, ax2 = plt.subplots()
    ax3 = ax2.twinx()
    ax2.plot(radius_set, velocity_radius_list, '--ok')
    ax2.grid()
    ax2.set_xlabel('r/R [m]')
    ax2.set_ylabel('Velocity [m/s]')
    fig2.savefig(f'{saveFigDir}/Velocity_R{i}.png')
    fig2.savefig(f'{saveFigDir}/Velocity_R{i}.pdf')
    plt.close(fig2)'''
    pbar.update(1)

v_t = pd.DataFrame(fluid.time_list)
v_t = pd.concat([v_t, pd.DataFrame(normalized_velocity_list)], axis = 1)
v_t.columns = ["time", "Velocity_ratio"]
v_t.to_csv(f"{currentPath}/Velocity_t.csv")

#Plot 
fig1, ax1 = plt.subplots()
ax1.plot(fluid.time_list, normalized_velocity_list, '-k')
#ax1.plot(fluid.time_list, np.repeat(1, len(fluid.time_list)), '--r')#, alpha = 0.2
#ax0.spines['right'].set_color('red')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel(r'Velocity ratio between center walls $\left [\frac{m/s}{m/s} \right ]$')
fig1.tight_layout()
plt.grid()
fig1.savefig(f'{saveFigDir}/Velocity_t.pdf')
fig1.savefig(f'{saveFigDir}/Velocity_t.png')
ax0 = ax1.twinx()
ax0.plot(fluid.time_list, normalized_voidfraction_list, '.-r')
ax0.set_ylabel(r'Velocity ratio between center walls', color = 'red')
ax0.yaxis.set_label_position("right")
ax0.yaxis.tick_right()
fig1.tight_layout()
fig1.savefig(f'{saveFigDir}/Velocity_t_eps.pdf')
fig1.savefig(f'{saveFigDir}/Velocity_t_eps.png')

#Create an empty list to store the velocity as a radial function
velocity_radius_list = []

#Create a set of radius
radius_set = list(set(slice['R']))


for r in radius_set:
    velocity_radius_list.append(np.mean(slice['velocity_X'][slice['R'] == r]))

#Store velocity as a radial function
fig2, ax2 = plt.subplots()
ax2.plot(radius_set, velocity_radius_list, 'ok')
ax2.set_xlabel('Radius [m]')
ax2.set_ylabel('Velocity [m/s]')
fig2.savefig(f'{saveFigDir}/Velocity_R{i}.png')
fig2.savefig(f'{saveFigDir}/Velocity_R{i}.pdf')




