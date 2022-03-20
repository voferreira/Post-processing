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

clip_results = False
smooth = True

if smooth:
    print('Smoothing, so not clipping')
    clip_results = False
#############################################################################

#############################################################################
'''Importing Libraries'''
from math import e, pi
from tqdm import tqdm
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

particle = 'Alumina'
scheme = 'Q1-Q1'

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
def basic_interpolation(now, nxt, frac):
    # generate a linear interpolation between now and the next step
    # must sort based on IDs
    now_id = np.argsort(now['ID'])
    nxt_id = np.argsort(nxt['ID'])
    interp = now.copy()

    vel = now['Velocity_Z'][now_id]*(1 - frac) + nxt['Velocity_Z'][nxt_id]*frac
    # breakpoint()
    interp.points = now.points[now_id]*(1 - frac) + nxt.points[nxt_id]*frac
    interp['Velocity_Z'] = vel

    # finally, generate particles
    ds_glyph_interp = interp.glyph(scale='Diameter', geom=sphere)
    ds_glyph_interp.active_scalars_name = 'Velocity_Z'
    return ds_glyph_interp

def create_particles(df_name):
    df_ = df_name

    # create spheres scaled to the diameter of the array "Diameter"
    df_glyph = df_.glyph(scale='Diameter', geom=sphere)

    # compute normalized velocity
    #df_glyph['Velocity_Norm'] = np.linalg.norm(df_glyph['Velocity'], axis=1)
    return df_glyph
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
bounds = pv.Cube(center = [0, -0.025, 0], x_length = 1.1, y_length = 0.05, z_length = 0.1)
pbar = tqdm(total=len(list_vtu), desc="Reading VTU")
for i in range(0, len(list_vtu)):
    #Read DF from VTU files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')
    exec(f'df_{i}[\'Velocity_Z\'] = df_{i}[\'Velocity\'][:, 0]')
    if clip_results:
        exec(f'df_{i} = df_{i}.clip_box(bounds, invert=False)')
    pbar.update()


#Select a data to apply the slice   
exec(f'df = df_{len(list_vtu)-2}')

#create a sphere with diameter 1
sphere = pv.Sphere(theta_resolution=10, phi_resolution=10)

#Choose the white background for the figure
pv.set_plot_theme("document")

#Create a plotter
plotter = pv.Plotter(off_screen=None)

plotter.camera.position = (0, 1.2, 0)
plotter.camera.focal_point = (-0.15, 0, 0)
plotter.camera.roll -= 90

#plotter.camera_position = 'xy'
#plotter.camera.roll += 90

df_glyph = create_particles(df)

#Add data to the plotter
plotter.add_mesh(df_glyph, cmap = 'turbo', scalars = 'Velocity_Z', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': True, 'position_x' : 0.60, 'position_y' :0.05})

plotter.open_movie(f"{saveFigDir}/particles_interpolate.mp4", framerate=40, quality=8) #20

n_substeps = 10
n_total = (len(time_list) - 1)*n_substeps
pbar = tqdm(total=n_total, desc="Writing Frames")
for i in range(1, len(list_vtu)-1):
    
    if smooth:
        df_prev = df
        exec(f'df = df_{i}')
        for j in range(n_substeps):
            interp = basic_interpolation(df_prev, df, j/n_substeps)
            plotter.mesh.overwrite(interp)
            plotter.add_text(f"Time: {time_list[i]}", name='Time_label')
            plotter.write_frame()
            pbar.update()
    else:
        exec(f'df = df_{i}')
        df_glyph = create_particles(df)
        df_glyph.active_scalars_name = 'Velocity_Z'
        plotter.mesh.overwrite(df_glyph)
        plotter.add_text(f"Time: {time_list[i]}", name='Time_label')
        plotter.write_frame()
        pbar.update()

pv.close()



