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

try:
    column_name_complement = sys.argv[3]

except:
    pass

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
inlet_velocity = fluid.prm_dict['w']

try:
    eps_comparison = pd.read_excel(f'D:/results/{particle_name}/Abstract_Results.xlsx')
except:
    eps_comparison = pd.read_excel(f'/mnt/d/results/{particle_name}/Abstract_Results.xlsx')
eps_exp = eps_comparison['Experimental'][eps_comparison['U'] == inlet_velocity].values
eps_RZ = eps_comparison['R-Z'][eps_comparison['U'] == inlet_velocity].values

g = abs(particle.prm_dict['gz'])       #m/s^2
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
eps_list = []
voidfraction_list = []
total_deltaP = []

#Create a list with all x values for the pressure takes
z_pressure_takes = np.arange(0.1, Hb + 0.01, 0.01)

deltaP_analytical = Analytical(Np, rp, rhol, rhop, g, Area)

if os.path.isdir(currentPath + '/P_z') == False:
    os.mkdir(currentPath + '/P_z')

pbar = tqdm(total = len(fluid.time_list)-1, desc="Processing data")
for i in range(len(fluid.time_list)):
    #Define "df" as last time step df
    exec(f'df = fluid.df_{i}') 
    
    #Adjust data format
    '''df = df.point_data_to_cell_data()'''

    #Find the voidfraction spatial average for the cells with voidfraction below 1
    df_voidfraction = pd.DataFrame(df['void_fraction'])
    voidfraction_ave = df_voidfraction[df_voidfraction[0] < 1].mean()
    voidfraction_list.append(voidfraction_ave[0])

    #Take the first slice of the domain at pos0
    pos0 = [0, 0, 0.1]
    slice0 = df.slice(normal=[0, 0, 1], origin = pos0)
    p0 = np.mean(slice0['pressure'])

    #Create empty lists to fill with x values and pressure as function of x values
    p_z = []
    z_list = []

    #Loop through z values above z0
    for j in range(len(z_pressure_takes)):
        z_list.append(z_pressure_takes[j]-pos0[2])
        slice_z = df.slice(normal=[0, 0, 1], origin = [0, 0, z_pressure_takes[j]])
        p_z.append((p0 - np.mean(slice_z['pressure']))*rhol)

    #Store the total pressure drop
    total_deltaP.append(p_z[-1])
    
    #Apply least-squares to find the porosity
    regr = linear_model.LinearRegression() #fit_intercept = 0
    z_list = np.array(z_list).reshape(-1, 1)
    p_z = np.array(p_z).reshape(-1, 1)
    model = regr.fit(z_list, p_z)
    r_sq = model.score(z_list, p_z)

    #Find linear portion of the graph and adjust a linear function to it
    x = z_list
    y = p_z

    """for j in range(len(y)-1, -1, -1):
        if y[j] >= deltaP_analytical*0.95:
                x = np.delete(x, j).reshape(-1, 1)
                y = np.delete(y, j).reshape(-1, 1)"""
    
    for j in range(len(y)-1):
        if y[j] >= deltaP_analytical*0.8:
                x = x[:j]
                y = y[:j]
                break
    
    model = regr.fit(x, y)

    #Calculate bed voidage by the slope of the pressure drop curve
    eps = 1-(model.coef_[0][0]/((rhop-rhol)*g))
    eps_list.append(eps)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    if boundary_condition == 'partial slip':
        fig1.suptitle(f'{particle_name}: {scheme}, time = {fluid.time_list[i]}-{boundary_condition}-Layer: {fluid.prm_dict["boundary layer thickness"]} m')
    else:
        fig1.suptitle(f'{particle_name}: {scheme}, time = {fluid.time_list[i]}')
    ax1.plot(z_list, p_z, 'ok', ms = 5, label ='Simulation')
    ax1.plot(z_list, np.repeat(deltaP_analytical, len(z_list)), '.r', ms = 2, label = 'Analytical')
    ax1.plot(x, y, '.g')
    ax1.plot(x, model.predict(x), '-b')
    #plt.text(17.5, 0, r'$\varepsilon =  (-dp/dz)$')
    ax1.grid()
    ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
    ax1.set_xlabel(r'$Height \/\ [m]$')
    ax1.set_xlim(0, 1.02)
    ax1.set_ylim(0, deltaP_analytical*1.30)
    ax1.legend()
    ax1.annotate(r'$\varepsilon (-dp/dz) = {:1.2}$'.format(eps), (x[round(len(x)/2)], y[round(len(y)/2)]), xytext=(0.65, 0.4), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.04), fontsize=14, horizontalalignment='right', verticalalignment='top')
    fig1.savefig(f'{saveFigDir}/P_z/P_z-{i}.png')
    fig1.savefig(f'{saveFigDir}/P_z/P_z-{i}.pdf')
    plt.close(fig1)

    pbar.update(1)

#Export the total pressure results as a function of time
csv = pd.DataFrame([fluid.time_list, total_deltaP], index=['time', 'deltaP']).transpose()
csv.to_csv(f'{saveFigDir}/deltaP_t.csv')

csv2 = pd.DataFrame([fluid.time_list, eps_list], index=['time', 'eps']).transpose()
csv2.to_csv(f'{saveFigDir}/eps_t.csv')

#Plot pressure vs z
fig0, ax0 = plt.subplots(1, 1)
if boundary_condition == 'partial slip':
    fig0.suptitle(f'{particle_name}: {scheme} - {boundary_condition}: Layer: {fluid.prm_dict["boundary layer thickness"]} m')
else:
    fig0.suptitle(f'{particle_name}: {scheme}')
ax0.plot(fluid.time_list, eps_list, 'ok', label = 'Simulation - pressure drop', ms = 5)
ax0.plot(fluid.time_list, voidfraction_list, 'og', label = 'Simulation - Average among voidage < 1 cells')
ax0.plot(fluid.time_list, np.repeat(eps_exp, len(fluid.time_list)), '--r', label = 'Experimental (+/- 5%)')
ax0.fill_between(fluid.time_list, np.repeat(eps_exp, len(fluid.time_list))*0.95, np.repeat(eps_exp, len(fluid.time_list))*1.05, color = 'r', alpha = 0.3)
ax0.plot(fluid.time_list, np.repeat(eps_RZ, len(fluid.time_list)), '.b', label = 'Richardson & Zaki')
ax0.legend()
ax0.grid()
#ax0.set_xlim(0, 35)
ax0.set_ylim(0, 1)
ax0.set_ylabel(r'$Bed \/\ voidage \/\ [-]$')
ax0.set_xlabel(r'$Time \/\ [s]$')
fig0.savefig(f'{saveFigDir}/eps_t.png')
fig0.savefig(f'{saveFigDir}/eps_t.pdf')


eps_ave = np.average(eps_list[-10:])
eps_std = np.std(eps_list[-10:])

print('Average among the last 10 seconds: ' + str(eps_ave))
print('Standard Deviation among the last 10 seconds: ' + str(eps_std))
print('Average delta P = {}')


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
            column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-{drag_model}-Lift'
        except:
            column_excel = f'{model_type}-Q{velocity_order}-Q{pressure_order}-{boundary_condition}-{drag_model}'
    try:
        column_excel = column_excel + f"-{fluid.prm_dict['grid type']}"

    except:
        pass

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


height = Vp*Np/ (Area * (1 - eps_ave))
height_exp = Vp*Np/ (Area * (1 - eps_exp[0]))

print(f"Bed height = {height} -> {abs(height_exp - height)/height_exp * 100} % deviation of height")

height_of_cells = 1.1 /(132)

number_of_cells_missing = abs(height_exp - height)/height_of_cells

print(f"Number of cells missing = {number_of_cells_missing}")
