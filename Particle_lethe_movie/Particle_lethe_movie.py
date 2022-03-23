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
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
#from sklearn import linear_model
#from scipy.interpolate import make_interp_spline
from Lethe_pyvista_tools import *

#import glob

#import os
import sys
#############################################################################

#############################################################################
'''Simulation properties'''

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath

particle_name = sys.argv[2]

particle = Lethe_pyvista_tools(currentPath, prm_file_name='liquid_fluidized_bed.prm')
particle.read_lethe_to_pyvista('result_particles.pvd', last = 201)

model_type = particle.prm_dict['vans model'].replace('model', '')
velocity_order = str(int(particle.prm_dict['velocity order']))
pressure_order = str(int(particle.prm_dict['pressure order']))

scheme = f'Model type {model_type} - Q{velocity_order}-Q{pressure_order}'

dp = particle.prm_dict['diameter']
rp = dp/2
rhop = particle.prm_dict['density particles']
Np = particle.prm_dict['number']
inlet_velocity = particle.prm_dict['u']

eps_comparison = pd.read_excel('D:/results/alginate/Abstract_Results.xlsx')
eps_exp = eps_comparison['Experimental'][eps_comparison['U'] == inlet_velocity].values
eps_RZ = eps_comparison['R-Z'][eps_comparison['U'] == inlet_velocity].values

g = abs(particle.prm_dict['gx'])       #m/s^2
Vnp = Np * 4/3 * pi * rp**3

rhol = particle.prm_dict['density']  #kg/m^3
nu = particle.prm_dict['kinematic viscosity'] #m^2/s
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

    vel = now['Velocity_[m/s]'][now_id]*(1 - frac) + nxt['Velocity_[m/s]'][nxt_id]*frac
    # breakpoint()
    interp.points = now.points[now_id]*(1 - frac) + nxt.points[nxt_id]*frac
    interp['Velocity_[m/s]'] = vel

    # finally, generate particles
    ds_glyph_interp = interp.glyph(scale='Diameter', geom=sphere)
    ds_glyph_interp.active_scalars_name = 'Velocity_[m/s]'
    return ds_glyph_interp

def create_particles(df_name):
    df_ = df_name

    # create spheres scaled to the diameter of the array "Diameter"
    df_glyph = df_.glyph(scale='Diameter', geom=sphere)

    # compute normalized velocity
    #df_glyph['Velocity_Norm'] = np.linalg.norm(df_glyph['Velocity'], axis=1)
    return df_glyph
#############################################################################

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

bounds = pv.Cube(center = [0, -0.025, 0], x_length = 1.1, y_length = 0.05, z_length = 0.1)
for i in range(0, len(particle.list_vtu)):
    #Read DF from VTU files
    exec(f'particle.df_{i}[\'Velocity_[m/s]\'] = particle.df_{i}[\'Velocity\'][:, 0]')
    if clip_results:
        exec(f'particle.df_{i} = particle.df_{i}.clip_box(bounds, invert=False)')


exec(f'df_glyph = create_particles(particle.df_{len(particle.list_vtu)-2})')

#Add data to the plotter
plotter.add_mesh(df_glyph, cmap = 'turbo', scalars = 'Velocity_[m/s]', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': True, 'position_x' : 0.60, 'position_y' :0.05})

plotter.open_movie(f"{saveFigDir}/particles_interpolate.mp4", framerate=40, quality=8) #20

n_substeps = 10
n_total = (len(particle.time_list) - 1)*n_substeps
pbar = tqdm(total=n_total, desc="Writing Frames")
df = particle.df_0
for i in range(1, len(particle.list_vtu)-1):
    if smooth:
        df_prev = df
        exec(f'df = particle.df_{i}')
        for j in range(n_substeps):
            interp = basic_interpolation(df_prev, df, j/n_substeps)
            plotter.mesh.overwrite(interp)
            plotter.add_text(f"Time: {particle.time_list[i]}", name='Time_label')
            plotter.write_frame()
            pbar.update()
    else:
        exec(f'df = particle.df_{i}')
        df_glyph = create_particles(df)
        df_glyph.active_scalars_name = 'Velocity_[m/s]'
        plotter.mesh.overwrite(df_glyph)
        plotter.add_text(f"Time: {particle.time_list[i]}", name='Time_label')
        plotter.write_frame()
        pbar.update()

pv.close()



