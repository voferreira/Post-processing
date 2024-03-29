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
saveFigDir = currentPath

try:
    particle_name = sys.argv[2]
except:
    if "Alginate" in currentPath or "alginate" in currentPath:
        particle_name = "Alginate"
    elif "Alumina" in currentPath or "alumina" in currentPath:
        particle_name = "Alumina"

fluid = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')
fluid.read_lethe_to_pyvista('result_.pvd', first = 1)
#fluid.read_lethe_to_pyvista_parallel('result_.pvd', last = 200, processors=2)

particle = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')

model_type = fluid.prm_dict['vans model'].replace('model', '')
velocity_order = str(int(fluid.prm_dict['velocity order']))
pressure_order = str(int(fluid.prm_dict['pressure order']))

drag_model = fluid.prm_dict['drag model']

boundary_condition = fluid.prm_dict['bc 0']

scheme = f'Model type {model_type} - Q{velocity_order}-Q{pressure_order}'

dp = particle.prm_dict['diameter']
rp = dp/2
Vp = 4/3*pi*rp**3
rhop = particle.prm_dict['density particles']
Np = particle.prm_dict['number']
inlet_velocity = fluid.prm_dict['u']

try:
    press_comparison = pd.read_excel(f'D:/results/{particle_name}/Abstract_Pressure.xlsx')
except:
    press_comparison = pd.read_excel(f'/mnt/d/results/{particle_name}/Abstract_Pressure.xlsx')


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
'''Functions'''
def Analytical(Np, rp, rhol, rhop, g, Area):
    Mp = Np*4/3*pi*rp**3*rhop
    delta_p = Mp*(rhop - rhol)*g/(rhop * Area)
    return delta_p
#############################################################################

#Define eps_list and voidfraction_list to append value for each time-step
total_deltaP = []

#Create a list with all x values for the pressure takes
x_pressure_takes = np.arange(-0.45, Hb/2 + 0.01, 0.01)

deltaP_analytical = Analytical(Np, rp, rhol, rhop, g, Area)

if os.path.isdir(currentPath + '/P_x') == False:
    os.mkdir(currentPath + '/P_x')

pbar = tqdm(total = len(fluid.time_list)-1, desc="Processing data")
for i in range(len(fluid.time_list)):
    #Define "df" as last time step df
    exec(f'df = fluid.df_{i}') 

    #Take the first slice of the domain at pos0
    pos0 = [-0.45, 0, 0]
    slice0 = df.slice(normal=[1, 0, 0], origin = pos0)
    p0 = slice0.integrate_data()['pressure']/slice0.area

    #Create empty lists to fill with x values and pressure as function of x values
    p_x = []
    x_list = []

    #Loop through z values above z0
    for j in range(len(x_pressure_takes)):
        x_list.append(x_pressure_takes[j]-pos0[0])
        slice_x = df.slice(normal=[1, 0, 0], origin = [x_pressure_takes[j], 0, 0])
        p_x.append((p0 - slice_x.integrate_data()['pressure']/slice_x.area)[0]*rhol)

    #Store the total pressure drop
    total_deltaP.append(p_x[-1])

    pbar.update(1)

#Export the total pressure results as a function of time
csv = pd.DataFrame([fluid.time_list, total_deltaP], index=['time', 'deltaP']).transpose()
csv.to_csv(f'{saveFigDir}/deltaP_t.csv')

print(f'Average delta P = {np.mean(total_deltaP)} Pa')

if write_to_excel:
    try:
        excel_file_path = f'D:/results/{particle_name}/Abstract_Pressure.xlsx'
        excel_pd_average = pd.read_excel(excel_file_path, sheet_name = 'Average', index_col= 0)
    except:
        excel_file_path = f'/mnt/d/results/{particle_name}/Abstract_Pressure.xlsx'
        excel_pd_average = pd.read_excel(excel_file_path, sheet_name = 'Average', index_col= 0)
    excel_pd_std = pd.read_excel(excel_file_path, sheet_name = 'Std', index_col= 0)

    excel_pd_std.drop(excel_pd_std.filter(regex="Unnamed"),axis=1, inplace=True)

    if boundary_condition == 'partial slip':
        column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-Layer: {fluid.prm_dict["boundary layer thickness"]} m'
    else:    
        try:
            lift = fluid.prm_dict["saffman lift force"]
            column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-{drag_model}-Lift'
        except:
            column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-{drag_model}'
    if column_excel not in excel_pd_average.columns:
        excel_pd_average[column_excel] = ''
    if column_excel not in excel_pd_std.columns:
        excel_pd_std[column_excel] = ''

    excel_pd_average[column_excel][excel_pd_average['U'] == fluid.prm_dict['u']] = np.mean(total_deltaP)
    excel_pd_std[column_excel][excel_pd_average['U'] == fluid.prm_dict['u']] = np.std(total_deltaP)
    
    with pd.ExcelWriter(excel_file_path, engine='openpyxl') as writer:
        excel_pd_average.to_excel(writer, sheet_name="Average")
        excel_pd_std.to_excel(writer, sheet_name = 'Std')  

