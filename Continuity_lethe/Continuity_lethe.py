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
write_to_excel = True
#############################################################################

#############################################################################
'''Importing Libraries'''
from math import e, pi
from re import U
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from tqdm import tqdm
from Lethe_pyvista_tools import *
from openpyxl import load_workbook

import os
import sys
#############################################################################

#############################################################################
'''Simulation properties'''

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]
saveFigDir = currentPath

all_paths = []

[ all_paths.append(os.path.join(currentPath, name)) for name in os.listdir(currentPath) if os.path.isdir(os.path.join(currentPath, name)) ]

try:
    particle_name = sys.argv[2]
except:
    if "Alginate" in currentPath or "alginate" in currentPath:
        particle_name = "Alginate"
    elif "Alumina" in currentPath or "alumina" in currentPath:
        particle_name = "Alumina"

#############################################################################

#############################################################################

NMFR_ave = []
NMFR_std = []

for k in range(len(all_paths)):

    if "output" in all_paths[k]:
        break

    fluid = Lethe_pyvista_tools(all_paths[k], prm_file_name='liquid_fluidized_bed.prm')
    fluid.read_lethe_to_pyvista('result_.pvd', first = 1)

    NMFR = []

    pbar = tqdm(total = len(fluid.time_list)-1, desc="Processing data")
    for i in range(len(fluid.time_list)):
        #Define "df" as last time step df
        exec(f'df = fluid.df_{i}') 
        
        #Adjust data format
        df = df.point_data_to_cell_data()

        #Take the first slice of the domain at pos0
        slice0 = df.slice(normal=[1, 0, 0], origin = [-0.55, 0, 0])
        slicef = df.slice(normal=[1, 0, 0], origin = [0.55, 0, 0])

        slice0_cell_info = slice0.compute_cell_sizes()
        slicef_cell_info = slicef.compute_cell_sizes()

        U0 = slice0['velocity'][:, 0]
        Uf = slicef['velocity'][:, 0]

        NMFR_i = (np.dot(U0, slice0_cell_info['Area']) - np.dot(Uf, slicef_cell_info['Area']))/np.dot(U0, slice0_cell_info['Area'])

        NMFR.append(NMFR_i)

        pbar.update(1)

    #Export the total pressure results as a function of time
    csv = pd.DataFrame([fluid.time_list, NMFR], index=['time', 'NMFR']).transpose()
    csv.to_csv(f'{all_paths[k]}/NMFR_t.csv')

    NMFR_ave.append(np.mean(NMFR))
    NMFR_std.append(np.std(NMFR))


csv2 = pd.DataFrame([NMFR_ave, NMFR_std], index=['NMFR_ave', 'NMFR_std']).transpose()
csv2.to_csv(f'{currentPath}/NMFR_ave.csv')
print("done!")
