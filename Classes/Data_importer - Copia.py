#############################################################################
"""Class with data importers"""
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
#Define class:

class Data_importer:
    def lethe_fluid(self, path):
        #Read name of files in .pvd file
        files = pd.read_csv(f'{path}/result_.pvd',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
        files = files.dropna()
        self.time_list = files['time'].tolist()
        list_vtu = files['vtu'].tolist()
        list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in list_vtu]

        #Read VTU data
        for i in range(len(list_vtu)):
            #Read DF from VTU files
            exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')

    
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

    def openfoam(self, path):
        #Take case path as argument and store at currentPath
        currentPath = sys.argv[1]

        #Define list of VTK files and time list:
        list_vtk = os.listdir(currentPath)

        try:
            list_vtk.remove('walls')
            list_vtk.remove('inlet')
            list_vtk.remove('outlet')
            list_vtk.remove('VTK')
        except:
            print('')

        #Create a list of time_steps
        time_list = [i.replace('.vtk', '').replace('CFD_', '') for i in list_vtk]
        time_list = [float(i) for i in time_list]
        time_list = [round(2*2.5*i*10e-5, 2) for i in time_list]

        #Sort list
        time_list, list_vtk = (list(t) for t in zip(*sorted(zip(time_list, list_vtk))))

        #Read VTK data
        for i in range(0, len(list_vtk)):
            #Read DF from VTK files
            exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtk[i]}\')')