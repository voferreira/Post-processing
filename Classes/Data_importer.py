#############################################################################
"""Class with data importers"""
#############################################################################

#############################################################################

#Import modules
import sys
import os
import pandas as pd
import pyvista as pv

#Define class:
class Data_importer():

    def __init__(self):

        #Take case path as argument and store at currentPath
        self.path = sys.argv[1]
        
    def lethe_save_path(self):
        self.path_save = self.path.replace('/output', '')
    
    def lethe_fluid(self):
        #call function
        self.lethe_save_path()
        
        #Read name of files in .pvd file
        files = pd.read_csv(f'{self.path}/result_.pvd',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
        files = files.dropna()
        self.time_list = files['time'].tolist()
        list_vtu = files['vtu'].tolist()
        list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in list_vtu]

        #Read VTU data
        for i in range(len(list_vtu)):
            #Read DF from VTU files
            exec(f'self.df_{i} = pv.read(f\'{self.path}/{list_vtu[i]}\')')

    def lethe_particle(self):
        #define saving path
        self.lethe_save_path()

        #Read name of files in .pvd file
        files = pd.read_csv(f'{self.path}/result_particles.pvd',sep='"',skiprows=6, usecols=[1, 5], names = ['time', 'vtu'])
        files = files.dropna()
        time_list = files['time'].tolist()
        list_vtu = files['vtu'].tolist()
        list_vtu = [i.replace('.pvtu', '.0000.vtu') for i in list_vtu]

        #Read VTU data
        for i in range(0, len(list_vtu)):
            #Read DF from VTU files
            exec(f'self.df_{i} = pv.read(f\'{self.path}/{list_vtu[i]}\')')
        

    def cfdem_fluid(self):
        #define save path
        self.path_save = self.path.replace('/VTK', '')

        #Define list of VTK files and time list:
        list_vtk = os.listdir(self.path)
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
            exec(f'self.df_{i} = pv.read(f\'{self.path}/{list_vtk[i]}\')')