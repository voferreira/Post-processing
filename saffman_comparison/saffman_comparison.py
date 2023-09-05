# importing libraries
import os
import openpyxl
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from lethe_pyvista_tools import *

test = False
test_case_n = 6
write_vtu_files = True


def listdirs(rootdir):
    for file in os.scandir(rootdir):
        if file.is_dir():
            if "output" in file.path:
                subdirec.append(file.path.replace("output", ""))
                print(f"FOUND SUBDIRECTORY: {file.path.replace('output', '')}")
            
            listdirs(file.path)
            
    return subdirec



def modify(case_path):

    print("\n")
    print(f"WORKING ON: {i}")
    print("\n")

    particles = lethe_pyvista_tools(case_path = case_path, prm_file_name = "liquid_fluidized_bed.prm")

    try:
        particles.read_lethe_to_pyvista(pvd_name = "result__particles.pvd")
    except:
        particles.read_lethe_to_pyvista(pvd_name = "result_particles.pvd")

    force_radius = []
    force_radius_std = []

    for j in range(len(particles.df)):
        dif_components = particles.df[j]["FemForce"][:, 1] + particles.df[j]["FemForce"][:, 2] + 10**(-12)

        unit_vector = np.divide((dif_components), (np.absolute(dif_components)))

        particles.df[j]["force_radius"] = (particles.df[j]["FemForce"][:, 1]**2 + particles.df[j]["FemForce"][:, 2]**2)**(1/2) * unit_vector

        force_radius.append(np.mean(particles.df[j]["force_radius"]))
        force_radius_std.append(np.std(particles.df[j]["force_radius"]))

    df = pd.DataFrame({'time': particles.time_list, 'force_r': force_radius, 'force_r_std': force_radius_std})
    df.to_excel(f'{i}/force_radius.xlsx', index = False)

    if write_vtu_files:
        particles.write_vtu(prefix = "mod_")



subdirec = []
subdirectories = listdirs(os.getcwd())

fails = []

force_radius_list_mean = []

if test:
    i = subdirectories[test_case_n]
    modify(i)
else:
    for i in subdirectories:
        try:
            modify(i)
            print("\n")
            print(f"FINISHED for: {i}")
            print("\n")
        except Exception as e:
            print(e)
            print("\n")
            print(f"FAILED for : {i}")
            print("\n")

            fails.append(i)

print("\n")
print(f"LIST OF FAILS:")
print(fails)
print("\n")