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
import matplotlib.tri as tri
from matplotlib import cm
import pyvista as pv

import glob

import os
import sys
from pyvista.plotting.renderer import CameraPosition

from scipy.linalg.special_matrices import dft
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

#############################################################################

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath.replace('\output', '')

#Define list of VTK files and time list:
listVTK = os.listdir(currentPath)
listVTK = [x for x in listVTK if "particles" not in x ]
listVTK = [x for x in listVTK if "pvd" not in x ]
listVTK = [x for x in listVTK if "pvtu" not in x ]

#Reorder listVTK
listVTK = sorted(listVTK)

#Create a list of time_steps
time_list = {i.replace('.vtu', '').replace('result_.', '').replace('.', '') for i in listVTK}
time_list = {float(i) for i in time_list}
time_list = sorted({round(i/10000000, 2) for i in time_list})

#Read VTK data
for i in range(0, len(listVTK)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{listVTK[i]}\')')

#Select a data to apply the slice   
exec(f'df = df_{len(listVTK)-1}')

#Choose the white background for the figure
pv.set_plot_theme("document")

#Convert Cell data to point data, and then do it the opposite way
#to spread the cell data to the entire cell region
df = df.cell_data_to_point_data()
df = df.point_data_to_cell_data()

#Select a slice of the domain
single_slice = df.slice(normal=[0, 1, 0], origin = [0, 0, 0])

#Create a plotter
plotter = pv.Plotter(off_screen=None)

#Add data to the plotter
plotter.add_mesh(single_slice, scalars = single_slice['void_fraction'], cmap = 'turbo', preference = 'cells', show_edges=False)

#Add text to the plotter
plotter.add_text(f'Voidfraction profile \n at time = {2} s')

#Show and save the screenshot of the data
#cpos (camera position) -> cross-section
plotter.show(screenshot=f'{saveFigDir}/test.png', cpos="zx")
