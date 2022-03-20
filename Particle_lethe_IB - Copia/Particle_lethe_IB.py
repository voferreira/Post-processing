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
def basic_interpolation(field, now, nxt, frac):
    # generate a linear interpolation between now and the next step
    # must sort based on IDs
    now_id = np.argsort(now['ID'])
    nxt_id = np.argsort(nxt['ID'])
    interp = now.copy()

    vel = now[field][now_id]*(1 - frac) + nxt[field][nxt_id]*frac
    # breakpoint()
    interp.points = now.points[now_id]*(1 - frac) + nxt.points[nxt_id]*frac
    interp[field] = vel

    # finally, generate particles
    ds_glyph_interp = interp.glyph(scale='Diameter', geom=sphere)
    ds_glyph_interp.active_scalars_name = field
    return ds_glyph_interp

def basic_interpolation_fluid(field, now, nxt, frac):
    # generate a linear interpolation between now and the next step
    # must sort based on IDs
    now_id = np.argsort(now['ID'])
    nxt_id = np.argsort(nxt['ID'])
    interp = now.copy()

    vel = now[field][now_id]*(1 - frac) + nxt[field][nxt_id]*frac
    # breakpoint()
    interp.points = now.points[now_id]*(1 - frac) + nxt.points[nxt_id]*frac
    interp[field] = vel

    # finally, generate particles
    ds_glyph_interp = interp.slice(normal=[0, 1, 0], origin = [0, 0, 0])
    ds_glyph_interp.active_scalars_name = field
    return ds_glyph_interp

def create_particles(df_name):
    df_ = df_name

    # create spheres scaled to the diameter of the array "Diameter"
    df_glyph = df_.glyph(scale='Diameter', geom=sphere)

    # compute normalized velocity
    #df_glyph['Velocity_Norm'] = np.linalg.norm(df_glyph['Velocity'], axis=1)
    return df_glyph

def create_arrows(df_name):
    df_ = df_name

    # create spheres scaled to the diameter of the array "Diameter"
    df_glyph = df_.glyph(orient="velocity", scale = "velocity_mag", geom=arrows, factor=0.04, tolerance=0.04)

    # compute normalized velocity
    #df_glyph['Velocity_Norm'] = np.linalg.norm(df_glyph['Velocity'], axis=1)
    return df_glyph
#############################################################################

#particles = Lethe_pyvista_tools('D:/results/lethe/alginate/q1-q1', 'first_case.prm')
#particles.path_output('output')
#particles.read_lethe_to_pyvista('result_particles.pvd', last = 20)
fluid = Lethe_pyvista_tools('D:/results/lethe/alginate/q1-q1', 'first_case.prm')
fluid.path_output('output')
fluid.read_lethe_to_pyvista('result_.pvd')

#for i in range(0, len(particles.list_vtu)):
    #exec(f'particles.df_{i}[\'Velocity_x\'] = particles.df_{i}[\'Velocity\'][:, 0]')
#    exec(f'particles.df_{i}[\'Velocity_mag\'] = (particles.df_{i}[\'Velocity\'][:, 0]**2 + particles.df_{i}[\'Velocity\'][:, 1]**2 + particles.df_{i}[\'Velocity\'][:, 2]**2)**(1/2)')

for i in range(0, len(fluid.list_vtu)):
    exec(f'fluid.df_{i}[\'velocity_mag\'] = (fluid.df_{i}[\'velocity\'][:, 0]**2 + fluid.df_{i}[\'velocity\'][:, 1]**2 + fluid.df_{i}[\'velocity\'][:, 2]**2)**(1/2)')
    exec(f'fluid.df_{i}[\'ID\'] = np.arange(0, len(fluid.df_{i}[\'velocity_mag\']))')
    exec(f'fluid.df_{i}[\'velocity_x\'] = fluid.df_{i}[\'velocity\'][:, 0]')

#create a sphere with diameter 1
#sphere = pv.Sphere(theta_resolution=20, phi_resolution=20)


#create the arrows
#arrows = pv.Arrow(tip_resolution=30, shaft_resolution=30)

#Choose the white background for the figure
pv.set_plot_theme("document")

#Create a plotter
plotter = pv.Plotter(off_screen=None)

#plotter.camera_position = ([(-5.592937351020412, -4.137862113540812, 2.192970961001668),
# (-4.177361512303915, -0.1970836549756708, -0.12578112629332203),
# (-0.9551137293493976, 0.26392817505908883, -0.13453506018233582)])
#plotter.camera.focal_point = (0, 0, 0)
#plotter.camera.roll += 180

#plotter.camera_position = 'xy'
#plotter.camera.roll += 90

plotter.camera.position = (0, 2.3, 0)
plotter.camera.focal_point = (0, 0, 0)
plotter.camera.roll -= 90

#df_glyph = create_particles(particles.df_0)
#df_fluid_glyph = create_arrows(fluid.df_100.threshold(value = 1.0, scalars = 'velocity_mag'))
single_slice = fluid.df_20.slice(normal=[0, 1, 0], origin = [0, 0, 0])

#Add data to the plotter
#plotter.add_mesh(box, color="cyan", show_edges=True, opacity = 0.1, line_width = 0.2, smooth_shading=True)
#PARTICLES = plotter.add_mesh(df_glyph, cmap = 'turbo', scalars = 'Velocity_mag', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': True, 'position_x' : 0.90, 'position_y' :0.60})
FLUID = plotter.add_mesh(single_slice, scalars= 'velocity_x', cmap = 'turbo', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': True, 'position_x' : 0.60, 'position_y' :0.05})#clim = [0, 2.78]
#ARROWS = plotter.add_mesh(df_fluid_glyph, clim = [0, 2.78],scalars = 'velocity_mag', cmap = 'turbo', opacity = 0.4)
plotter.show_axes()

#plotter.open_gif(f"{particles.path_case}/particles_fluid2.gif") #20 #3.8
plotter.open_movie(f"{fluid.path_case}/fluid_vel.mp4", quality = 9, framerate = 40)

n_substeps = 10
n_total = (len(fluid.time_list) - 1)*n_substeps
pbar = tqdm(len(fluid.list_vtu)-1, desc="Writing Frames")
df = fluid.df_0
for i in range(2, len(fluid.list_vtu)):
    
    if smooth:
        exec(f'df_fluid = fluid.df_{i}')
        exec(f'df_prev = fluid.df_{i-1}')

        for j in range(n_substeps):
            interp = basic_interpolation_fluid('velocity_x', df_prev, df_fluid, j/n_substeps)
            plotter.mesh.overwrite(interp)
            plotter.add_text(f"Time: {fluid.time_list[i]} sec.", name='Time_label')
            plotter.write_frame()
        pbar.update()
    else:
        plotter.remove_actor(PARTICLES)
        plotter.remove_actor(FLUID)
        #exec(f'df = particles.df_{i}')
        exec(f'df_fluid = fluid.df_{i}')
        #df_glyph = create_particles(df)

        single_slice = df_fluid.slice(normal=[0, 1, 0], origin = [0, 0, 0])
        FLUID = plotter.add_mesh(single_slice, scalars= 'void_fraction', clim = [0.4, 1], cmap = 'turbo', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': True, 'position_x' : 0.90, 'position_y' :0.60})

        #df_glyph.active_scalars_name = 'Velocity_mag'

        #PARTICLES = plotter.add_mesh(df_glyph, cmap = 'turbo', scalars = 'Velocity_mag', smooth_shading=True, scalar_bar_args={'color': 'k', 'vertical': True}) #clim = [0, 2.78]
    #    plotter.mesh.overwrite(df_glyph)
        plotter.add_text(f"Time: {particles.time_list[i]} sec.", name='Time_label')
        plotter.write_frame()
        pbar.update()
plotter.close()


