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
#from sklearn import linear_model
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

try:
    column_name_complement = sys.argv[3]

except:
    pass

fluid = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')
#fluid.read_lethe_to_pyvista('result_.pvd', first = 1)
#fluid.read_lethe_to_pyvista_parallel('result_.pvd', last = 200, processors=2)

particle = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')
particle.read_lethe_to_pyvista('result__particles.pvd', first = 1, interval = 100)

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

X_acc_list = Y_acc_list = Z_acc_list = []


pbar = tqdm(total = len(particle.df_0['ID'])*(len(particle.time_list)-1), desc="Processing data")

for i in range(1, len(particle.time_list)):
    exec(f'df_0 = particle.df_{i - 1}')
    exec(f'df_1 = particle.df_{i}')
    
    X_list = []
    Y_list = []
    Z_list = []

    X_acc = 0
    Y_acc = 0
    Z_acc = 0

    for j in range(len(df_0['ID'])):
        pos_0 = df_0.points[j]
        pos_1 = df_1.points[df_1['ID'] == df_0['ID'][j]][0]

        signal_0 = np.sign(pos_0)
        signal_1 = np.sign(pos_1)

        if signal_0[0] == signal_1[0]:
            X_list.append(1)
            X_acc += 1
        else:
            X_list.append(0)

        if signal_0[1] == signal_1[1]:
            Y_list.append(1)
            Y_acc += 1
        else:
            Y_list.append(0)

        if signal_0[2] == signal_1[2]:
            Z_list.append(1)
            Z_acc += 1
        else:
            Z_list.append(0)

        pbar.update(1)
    print(X_acc)
    exec(f'particle.df_{i}["X_list"] = {X_list}')
    exec(f'particle.df_{i}["Y_list"] = {Y_list}')
    exec(f'particle.df_{i}["Z_list"] = {Z_list}')

    X_acc_list.append(X_acc)
    Y_acc_list.append(Y_acc)
    Z_acc_list.append(Z_acc)

print(X_acc_list)

plt.plot(particle.time_list, X_acc_list)
plt.show()

exit()

if write_to_excel:
    try:
        excel_file_path = f'D:/results/{particle_name}/Abstract_Results.xlsx'
        excel_pd_average = pd.read_excel(excel_file_path, sheet_name = 'Average', index_col= 0)
    except:
        excel_file_path = f'/mnt/d/results/{particle_name}/Abstract_Results.xlsx'
        excel_pd_average = pd.read_excel(excel_file_path, sheet_name = 'Average', index_col= 0)
    excel_pd_std = pd.read_excel(excel_file_path, sheet_name = 'Std', index_col= 0)

    excel_pd_std.drop(excel_pd_std.filter(regex="Unnamed"),axis=1, inplace=True)

    if boundary_condition == 'partial slip':
        column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-Layer: {fluid.prm_dict["boundary layer thickness"]} m'
    else:
        try:
            lift = fluid.prm_dict["saffman lift force"]
            column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-{drag_model}-Lift-{fluid.prm_dict["grid type"]}'
        except:
            column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-{drag_model}-{fluid.prm_dict["grid type"]}'
    if column_name_complement:
        column_excel = column_excel + f"-{column_name_complement}"
    if column_excel not in excel_pd_average.columns:
        excel_pd_average[column_excel] = ''
    if column_excel not in excel_pd_std.columns:
        excel_pd_std[column_excel] = ''

    excel_pd_average[column_excel][excel_pd_average['U'] == fluid.prm_dict['u']] = eps_ave
    excel_pd_std[column_excel][excel_pd_average['U'] == fluid.prm_dict['u']] = eps_std
    
    with pd.ExcelWriter(excel_file_path, engine='openpyxl') as writer:
        excel_pd_average.to_excel(writer, sheet_name="Average")
        excel_pd_std.to_excel(writer, sheet_name = 'Std')  



