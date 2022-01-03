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
scheme = 'Q2-Q1'

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
files = pd.read_csv(f'{currentPath}/result_particles.pvd',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
files = files.dropna()
time_list = files['time'].tolist()
list_vtu = files['vtu'].tolist()
list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in list_vtu]

#Read VTU data
for i in range(0, len(list_vtu)):
    #Read DF from VTU files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')


#Select a data to apply the slice   
exec(f'df = df_{len(list_vtu)-1}')
df = df
print(df.points[0])
#Choose the white background for the figure
pv.set_plot_theme("document")

#Create a plotter
plotter = pv.Plotter(off_screen=None)

#Apply spherical glyph to the particles:
spheres = pv.Sphere(radius = rp, center = df.points.tolist(), theta_resolution=30, phi_resolution=30)#df.glyph(scale = 'Velocity', factor = 1, progress_bar = True)

#Add data to the plotter
plotter.add_mesh(spheres, cmap = 'turbo')#, scalars = df['Velocity'], cmap = 'turbo', render_points_as_spheres = True, point_size = 5)#df['Diameter'][0])

plotter.camera_position = 'xy'
plotter.camera.roll += 90

plotter.show(screenshot=f'{saveFigDir}/particles.png')
