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

particle = 'Alginate'
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
saveFigDir = currentPath.replace('/output', '')

#Read name of files in .pvd file
files = pd.read_csv(f'{currentPath}/result_.pvd',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
files = files.dropna()
time_list = files['time'].tolist()
list_vtu = files['vtu'].tolist()
list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in list_vtu]

#Read VTU data
for i in range(0, len(list_vtu)):
    #Read DF from VTU files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')

#Define eps_list and voidfraction_list to append value for each time-step
eps_list = []
voidfraction_list = []
total_deltaP = []

#Create a list with all x values for the pressure takes
x_pressure_takes = np.arange(-0.45, Hb/2 + 0.01, 0.01)

deltaP_analytical = Analytical(Np, rp, rhol, rhop, g, Area)

for i in range(len(list_vtu)):
    #Define "df" as last time step df
    exec(f'df = df_{i}') 
    
    #Adjust data format
    df = df.point_data_to_cell_data()

    #Find the voidfraction spatial average for the cells with voidfraction below 1
    df_voidfraction = pd.DataFrame(df['void_fraction'])
    voidfraction_ave = df_voidfraction[df_voidfraction[0] < 1].mean()
    voidfraction_list.append(voidfraction_ave[0])

    #Take the first slice of the domain at pos0
    pos0 = [-0.45, 0, 0]
    slice0 = df.slice(normal=[1, 0, 0], origin = pos0)
    p0 = np.mean(slice0['pressure'])

    #Create empty lists to fill with x values and pressure as function of x values
    p_x = []
    x_list = []

    #Loop through z values above z0
    for j in range(len(x_pressure_takes)):
        x_list.append(x_pressure_takes[j]-pos0[0])
        slice_x = df.slice(normal=[1, 0, 0], origin = [x_pressure_takes[j], 0, 0])
        p_x.append((p0 - np.mean(slice_x['pressure']))*rhol)

    #Store the total pressure drop
    total_deltaP.append(-p_x[-1])
    
    #Apply least-squares to find the porosity
    regr = linear_model.LinearRegression()
    x_list = np.array(x_list).reshape(-1, 1)
    p_x = np.array(p_x).reshape(-1, 1)
    model = regr.fit(x_list, p_x)
    r_sq = model.score(x_list, p_x)

    #Find linear portion of the graph and adjust a linear function to it
    x = x_list
    y = p_x

    """for j in range(len(y)-1, -1, -1):
        if y[j] >= deltaP_analytical*0.95:
                x = np.delete(x, j).reshape(-1, 1)
                y = np.delete(y, j).reshape(-1, 1)"""
    
    for j in range(len(y)-1):
        if y[j] >= deltaP_analytical*0.95:
                x = x[:j]
                y = y[:j]
                break
    
    model = regr.fit(x, y)

    #Calculate bed voidage by the slope of the pressure drop curve
    eps = 1-(model.coef_[0][0]/((rhop-rhol)*g))
    eps_list.append(eps)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    fig1.suptitle(f'{particle}: {scheme}, time = {time_list[i]}')
    ax1.plot(x_list, p_x, 'ok', ms = 5)
    ax1.plot(x, y, '.g')
    ax1.plot(x, model.predict(x), '-b')
    #plt.text(17.5, 0, r'$\varepsilon =  (-dp/dz)$')
    ax1.grid()
    ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
    ax1.set_xlabel(r'$Height \/\ [m]$')
    ax1.set_xlim(0, 1.02)
    ax1.set_ylim(0, deltaP_analytical*1.30)
    ax1.annotate(r'$\varepsilon (-dp/dz) = {:1.2}$'.format(eps), (x[round(len(x)/2)], y[round(len(y)/2)]), xytext=(0.65, 0.4), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.04), fontsize=14, horizontalalignment='right', verticalalignment='top')
    fig1.savefig(f'{saveFigDir}/P_x/P_x-{i}.png')
    plt.close(fig1)

#Export the total pressure results as a function of time
csv = pd.DataFrame([time_list, total_deltaP], index=['time', 'deltaP']).transpose()
csv.to_csv(f'{saveFigDir}/deltaP_t.csv')

#Determine a moving average of the eps values
eps_ave = np.mean(eps_list[-22:])
eps_moving_ave = []
n_moving_ave = 4


#Plot pressure vs z
fig0, ax0 = plt.subplots(1, 1)
fig0.suptitle(f'{particle}: {scheme}')
ax0.plot(time_list, eps_list, 'ok', label = 'Simulation - pressure drop', ms = 5)
#ax0.plot(time_list, voidfraction_list, 'og', label = 'Simulation - Average among voidage < 1 cells')
ax0.plot(time_list, np.repeat(eps_exp, len(time_list)), '--r', label = 'Experimental (+/- 10%)')
ax0.fill_between(time_list, np.repeat(eps_exp, len(time_list))*0.95, np.repeat(eps_exp, len(time_list))*1.05, color = 'r', alpha = 0.3)
ax0.plot(time_list, np.repeat(eps_RZ, len(time_list)), '.b', label = 'Richardson & Zaki')
ax0.legend()
ax0.grid()
ax0.set_ylabel(r'$Bed \/\ voidage \/\ [-]$')
ax0.set_xlabel(r'$Time \/\ [s]$')
fig0.savefig(f'{saveFigDir}/eps_t.png')

# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# fig1.suptitle(f'{particle}: {scheme}')
# ax1.plot(x_list, p_x, 'ok', ms = 5)
# ax1.plot(x, y, '.g')
# ax1.plot(x, model.predict(x), '-b')
# #plt.text(17.5, 0, r'$\varepsilon =  (-dp/dz)$')
# ax1.grid()
# ax1.set_ylabel(r'$\Delta P \/\ [Pa]$')
# ax1.set_xlabel(r'$Height \/\ [m]$')
# fig1.savefig(f'{saveFigDir}/P_x.png')

height = Vnp/(Area*(1-eps_ave))

print('The average porosity of the bed is:\n')
print('eps = {}'.format(eps_ave), 'height = {}'.format(height))

exit()
