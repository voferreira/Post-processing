#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira
Date: 2021/11/03
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
import paramiko

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

#Access cluster
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    ssh.connect(hostname = '200.136.230.156', username = 'victor', password = 'victorcluster', port=22)
except:
    print('Connection failed!')
    print('Exiting...')
    exit()

print('Connection established')

#Establishing SFTP to access simulation data 
sftp_client = ssh.open_sftp()
sftp = paramiko.SFTPClient

#Access simulation directory
simulationName = '/home/victor/CFDEM/victor-PUBLIC-5.x/run/Artigo1/Alginate/new'
sftp_client.chdir(simulationName)

#Reconstructing OpenFOAM parallel case
"""stdin, stdout, stderr = ssh.exec_command(f"cd {simulationName}/CFD;reconstructPar -noLagrangian")
stdout = stdout.readlines()
stdout = "".join(stdout)
print(f'Output of running \'reconstructPar -noLagrangian\' in {simulationName}/CFD:')
print(stdout)"""


"""stdin, stdout, stderr = ssh.exec_command(f"source /home/victor/.bashrc;cd {simulationName}/CFD; echo $HOME; echo $CFDEM_VERSION", get_pty = True)
print(f'Output of running \'foamToVTK\' in {simulationName}/CFD:')
stdout = stdout.readlines()
stdout = "".join(stdout)
print(stdout)"""


#Find names of VTK files
stdin, stdout, stderr = ssh.exec_command(f"cd {simulationName}/CFD/VTK; ls -p | grep -v / ", get_pty = True)
print(f'Reading following files in {simulationName}:')
stdout = stdout.readlines()
stdout = "".join(stdout)
filesVTK = stdout.split()

filesVTK = {i.replace('.vtk', '').replace('CFD_', '') for i in filesVTK}

print(filesVTK)

#Ploting voidfraction as a func of time
for f in filesVTK:
    #Reading DF from VTK files
    exec(f'df_{f} = pv.read(f\'{currentPath}/CFD/VTK/CFD_{f}.vtk\')')
    
    #Converting voidfraction data to pandas df
    exec(f'df_{f}_voidfraction = pd.DataFrame(df_{f}.cell_data[\'voidfraction\'])')
    
    #printing voidfraction
    #exec(f'print(df_{time}_voidfraction[df_{time}_voidfraction != 1].mean()[0])')
    
    if figure:
        #ploting
        exec(f'plt.scatter({f}*0.0002, df_{f}_voidfraction[df_{f}_voidfraction != 1].mean()[0], c = \'k\')')

if figure:
    plt.show()

#Close connection to cluster
sftp_client.close()
ssh.close()

exit()

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