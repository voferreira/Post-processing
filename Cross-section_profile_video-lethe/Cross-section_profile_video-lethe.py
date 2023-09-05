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
import cv2

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

field_name = 'velocity'

#############################################################################

#############################################################################
'''Functions'''

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

#Read VTK data
for i in range(0, len(list_vtu)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')
    exec(f'df_{i}.velocity_x = df_{i}["velocity"][:, 0]')

#Selecet a data to apply the slice   
exec(f'df = df_{0}')

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

plotter.camera.position = (0, 2.5, 0.55)
plotter.camera.focal_point = (0, 0, 0.55)
plotter.camera.roll -= 90

#Create a movie
video_name = f'{saveFigDir}/{field_name}.mp4'
plotter.open_movie(video_name, framerate = 3.8, quality = 5)

#Add data to the plotter
scale_bar = dict(vertical=True)
#plotter.add_mesh(single_slice, scalars = single_slice['voidfraction'], cmap = 'turbo', preference = 'cells', show_edges=False) #flip_scalars=True

#Add text to the plotter
#plotter.add_text(f'Time = {5} s')

#Show and save the screenshot of the data
#cpos (camera position) -> cross-section
#plotter.show(screenshot=f'{saveFigDir}/test.png', cpos="xz")

# Run through each frame
#plotter.write_frame()  # write initial data

# Update scalars on each frame
for i in range(1, len(list_vtu)):
    exec(f'df = df_{i}')
    #df = df.cell_data_to_point_data()
    #df = df.point_data_to_cell_data()
    single_slice = df.slice(normal=[0, 1, 0], origin = [0, 0, 0])
    plotter.add_mesh(single_slice, scalars = single_slice[field_name], cmap = 'turbo', show_edges=False, scalar_bar_args = scale_bar) # preference = 'cells', flip_scalars=True
    plotter.write_frame()  # Write this frame

# Be sure to close the plotter when finished
plotter.close()
