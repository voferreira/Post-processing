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
smooth = False

if smooth:
    print('Smoothing, so not clipping')
    clip_results = False
#############################################################################

#############################################################################
'''Importing Libraries'''
from Lethe_pyvista_tools import *

from math import pi
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
from sklearn import linear_model
from scipy.interpolate import make_interp_spline

#############################################################################

#############################################################################
'''Functions'''
def basic_interpolation(now, nxt, frac):
    # generate a linear interpolation between now and the next step
    # must sort based on IDs
    now_id = np.argsort(now['ID'])
    nxt_id = np.argsort(nxt['ID'])
    interp = now.copy()

    vel = now['Velocity_x'][now_id]*(1 - frac) + nxt['Velocity_x'][nxt_id]*frac
    # breakpoint()
    interp.points = now.points[now_id]*(1 - frac) + nxt.points[nxt_id]*frac
    interp['Velocity_x'] = vel

    # finally, generate particles
    ds_glyph_interp = interp.glyph(scale='Diameter', geom=sphere)
    ds_glyph_interp.active_scalars_name = 'Velocity_x'
    return ds_glyph_interp

def create_particles(df_name):
    df_ = df_name

    # create spheres scaled to the diameter of the array "Diameter"
    df_glyph = df_.glyph(scale='Diameter', geom=sphere)

    # compute normalized velocity
    #df_glyph['Velocity_Norm'] = np.linalg.norm(df_glyph['Velocity'], axis=1)
    return df_glyph

#############################################################################

particles = Lethe_pyvista_tools('D:/multiple particles sedimentation', 'IB_flow_PBR.prm')
particles.path_output('')
particles.read_lethe_to_pyvista('ib_particles_reconstruct.pvd')

for i in range(0, len(particles.list_vtu)):
    exec(f'particles.df_{i}[\'Velocity_x\'] = particles.df_{i}[\'Velocity\'][:, 0]')
    


#create a sphere with diameter 1
sphere = pv.Sphere(theta_resolution=50, phi_resolution=50)

#create the box
box = pv.Box(bounds=(-5.25, -3.25, - 0.25, 0.25, - 0.25, 0.25))


#Choose the white background for the figure
pv.set_plot_theme("document")

#Create a plotter
plotter = pv.Plotter(off_screen=None)

plotter.camera_position = ([(-5.592937351020412, -4.137862113540812, 2.192970961001668),
 (-4.177361512303915, -0.1970836549756708, -0.12578112629332203),
 (-0.9551137293493976, 0.26392817505908883, -0.13453506018233582)])
#plotter.camera.focal_point = (0, 0, 0)
#plotter.camera.roll += 180

#plotter.camera_position = 'xy'
#plotter.camera.roll += 90

df_glyph = create_particles(particles.df_1000)
df_fluid_glyph = create_arrows()

#Add data to the plotter
plotter.add_mesh(box, color="cyan", show_edges=True, opacity = 0.1, line_width = 0.2, smooth_shading=True)
plotter.add_mesh(df_glyph, cmap = 'turbo', scalars = 'Velocity_x', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': False})
plotter.show_axes()
plotter.open_gif(f"{particles.path_case}/particles_fluid.gif") #20 #3.8

n_substeps = 5
n_total = (len(particles.time_list) - 1)*n_substeps
pbar = tqdm(len(particles.list_vtu), desc="Writing Frames")
df = particles.df_0
for i in range(1, len(particles.list_vtu), 10):
    
    if smooth:
        df_prev = df
        exec(f'df = particles.df_{i}')
        for j in range(n_substeps):
            interp = basic_interpolation(df_prev, df, j/n_substeps)
            plotter.mesh.overwrite(interp)
            plotter.add_text(f"Time: {particles.time_list[i]}", name='Time_label')
            plotter.write_frame()
        pbar.update()
    else:
        exec(f'df = particles.df_{i}')
        df_glyph = create_particles(df)
        df_glyph.active_scalars_name = 'Velocity_x'
        plotter.mesh.overwrite(df_glyph)
        plotter.add_text(f"Time: {particles.time_list[i]}", name='Time_label')
        plotter.write_frame()
        pbar.update()

pv.close()



