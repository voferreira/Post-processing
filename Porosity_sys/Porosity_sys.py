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
from math import e, pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from sklearn import linear_model
from scipy.interpolate import make_interp_spline

import glob

import os
import sys
#############################################################################

#############################################################################
'''Simulation properties'''

g = 9.81        #m/s^2
rp = 0.00154362 #            Alumina: 0.00154362 Alginate: 0.001332  #m
dp = 2*rp       #m
rhop = 3585.892  #         Alumina: 3585.892   Alginate: 1028.805 #kg/m^3
Np = 72401 #               Alumina: 72401      Alginate: 107960    #[-]
Vnp = Np * 4/3 * pi * rp**3

rhol = 996.778  #kg/m^3
nu = 0.00000084061 #m^2/s
mu = rhol*nu    #Pa
Db = 0.1        #m
Rb = Db/2       #m
Area = pi*Rb**2 #m^2

inlet_velocity = 0.0951 #m/s

eps_RZ = 0.550
eps_exp = 0.550

#############################################################################

#############################################################################
'''Functions'''
def Analytical(Np, rp, rhol, rhop, g, Area):
    Mp = Np*4/3*pi*rp**3*rhop
    delta_p = Mp*(rhop - rhol)*g/(rhop * Area)
    return delta_p
#############################################################################

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath.replace('/VTK', '')

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]

#Define list of VTK files and time list:
listVTK = os.listdir(currentPath)

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

#Read VTK data
for i in range(0, len(listVTK)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{listVTK[i]}\')')

#Define eps_list and voidfraction_list to append value for each time-step
eps_list = []
voidfraction_list = []

deltaP_analytical = Analytical(Np, rp, rhol, rhop, g, Area)

for i in range(len(listVTK)):
    #Define "df_last" as last time step df
    exec(f'df = df_{i}') 

    #Convert cell data to point data
    df = df.cell_centers()

    #Convert cell positions, pressure, and voidfraction values to pandas
    df_cells = pd.DataFrame(df.points)
    df_p     = pd.DataFrame(df['p'])
    df_voidfraction = pd.DataFrame(df['voidfraction'])

    #Find the voidfraction spatial average for the cells with voidfraction below 1
    voidfraction_ave = df_voidfraction[df_voidfraction[0] < 1].mean()
    voidfraction_list.append(voidfraction_ave[0])
    
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
    #model = regr.fit(z_list, p_z)
    #r_sq = model.score(z_list, p_z)

    #Find linear portion of the graph and adjust a linear function to it
    x = z_list
    y = p_z

    for j in range(len(y)-1, -1, -1):
        if y[j] >= deltaP_analytical*0.95:
                x = np.delete(x, j).reshape(-1, 1)
                y = np.delete(y, j).reshape(-1, 1)
    
    model = regr.fit(x, y)

    #Calculate bed voidage by the slope of the pressure drop curve
    eps = 1-(model.coef_[0][0]/((rhop-rhol)*g))
    eps_list.append(eps)

#Determine a moving average of the eps values
"""eps_ave = np.mean(eps_list[-22:])
eps_moving_ave = []
n_moving_ave = 4

for i in range(n_moving_ave, len(eps_list)):
    aux = 0
    for j in range(n_moving_ave):
        aux += eps_list[i-j]
    eps_moving_ave.append(aux/n_moving_ave)"""

#Smooth the moving average values:
#time_eps_spline = make_interp_spline(time_list[n_moving_ave:], eps_moving_ave)
#time_smooth_moving_ave = np.linspace(time_list[n_moving_ave], time_list[-1], 500)
#eps_smooth_moving_ave = time_eps_spline(time_smooth_moving_ave)

#Plot voidage as a function of time
fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
fig0.suptitle('Explicit scheme')
ax0.plot(time_list, eps_list, 'ok', label = 'Simulation - pressure drop', ms = 5)
ax0.plot(time_list, np.repeat(eps_exp, len(time_list)), '--r', label = 'Experimental (+/- 10%)')
ax0.fill_between(time_list, np.repeat(eps_exp, len(time_list))*0.95, np.repeat(eps_exp, len(time_list))*1.05, color = 'r', alpha = 0.3)
ax0.plot(time_list, np.repeat(eps_RZ, len(time_list)), '.b', label = 'Richardson & Zaki')
#plt.plot(time_list[n_moving_ave:], eps_moving_ave, '--b', label = f'Simulation moving average (last {n_moving_ave} points)')
ax0.plot(time_list, voidfraction_list, 'og', label = 'Simulation - Average among voidage < 1 cells')
#ax0.plot(time_smooth_moving_ave, eps_smooth_moving_ave, '--b', label = f'Simulation - Moving average (last {n_moving_ave} points)')
#plt.text(17.5, 0, r'$\varepsilon =  (-dp/dz)$')
ax0.legend()
ax0.grid()
ax0.set_ylabel(r'$Bed \/\ voidage \/\ [-]$')
ax0.set_xlabel(r'$Time \/\ [s]$')
fig0.savefig(f'{saveFigDir}/eps_t.png')

#plot p as a function o z
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig1.suptitle('Explicit scheme')
ax1.plot(z_list, p_z, 'ok', ms = 5)
ax1.plot(x, y, '.g')
ax1.plot(x, model.predict(x), '-b')
#plt.text(17.5, 0, r'$\varepsilon =  (-dp/dz)$')
ax1.grid()
ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
ax1.set_xlabel(r'$Height \/\ [m]$')
fig1.savefig(f'{saveFigDir}/P_z.png')


height = Vnp/(Area*(1-eps_ave))

print('The average porosity of the bed is:\n')
print('eps = {}'.format(eps_ave), 'height = {}'.format(height))

exit()
