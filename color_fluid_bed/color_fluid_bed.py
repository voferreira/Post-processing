# importing libraries
import os
from time import time

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from lethe_pyvista_tools import *

from lethe_pyvista_tools import *

test = False
first = 0
test_case_n = 2
write_vtu_files = False
drop_duplicates_true = False
n_procs = 2
folder_components = ["output"]


def listdirs(rootdir):
    for file in os.scandir(rootdir):
        if file.is_dir():
            present_in_file = [elem in file.path for elem in folder_components]
            if all(present_in_file):
                subdirec.append(file.path.replace("output", ""))
                print(f"FOUND SUBDIRECTORY: {file.path.replace('output', '')}")
            
            listdirs(file.path)
            
    return subdirec



def modify(case_path):

    if "alginate" in case_path:
        #first = 60
        reference_time_step = 20
        step = 2
        criterion = 0.5
        last = None
    elif "/e" in case_path:
        #first = 50
        reference_time_step = 25
        step = 2
        criterion = 0.95
        last = 202
        print("NEW \n")
    elif "alumina" in case_path:
        #first = 30
        step = 2
        reference_time_step = 25
        criterion = 0.95
        last = None

    print("\n")
    print(f"WORKING ON: {i}")
    print("\n")

    particles = lethe_pyvista_tools(case_path = case_path, prm_file_name = "liquid_fluidized_bed.prm", pvd_name = "result__particles.pvd", step = step, last = last, n_procs = 4)

    #particles.read_lethe_to_pyvista(pvd_name = "result__particles.pvd", first = first)

    #reference_time_step = 0#30

    split_height = np.mean(particles.get_df(reference_time_step).points[:, 0])
    split_radius = 0.05/(2**(1/2))

    particles.modify_array(array_name = "height_color", condition = f"x > {split_height}", array_values = 1, restart_array = False, reference_time_step = reference_time_step)

    particles.modify_array(array_name = "radius_color", condition = f"(z**2 + y**2)**(1/2) > {split_radius}", array_values = 1, restart_array = False, reference_time_step = reference_time_step)

    particles.modify_array(array_name = "theta_color", condition = "z > 0", array_values = 1, restart_array = False, reference_time_step = reference_time_step)

    n_ID_height = []
    n_ID_radius = []
    n_ID_theta = []

    n_flows_height = -1
    n_flows_radius = -1
    n_flows_theta = -1
    n_flows_doucet = -1

    time_height = -1
    time_radius = -1
    time_theta = -1
    time_doucet = -1

    flow_throughs = np.multiply(particles.prm_dict.get("u"), np.subtract(particles.time_list, particles.time_list[reference_time_step]))

    particles.mixing_index_nearest_neighbors(reference_array = "height_color", mixing_index_array_name = "height_mixing_index")
    height_mixing_index = particles.mixing_index
    height_mixing_index_std = particles.mixing_index_std
    particles.mixing_index_nearest_neighbors(reference_array = "radius_color", mixing_index_array_name = "radius_mixing_index")
    radius_mixing_index = particles.mixing_index
    radius_mixing_index_std = particles.mixing_index_std
    particles.mixing_index_nearest_neighbors(reference_array = "theta_color", mixing_index_array_name = "theta_mixing_index")
    theta_mixing_index = particles.mixing_index
    theta_mixing_index_std = particles.mixing_index_std
    particles.mixing_index_doucet(reference_time_step = reference_time_step, use_cyl = True, increasing_index = True, normalize = True)
    doucet_mixing_index = particles.mixing_index
    doucet_eigenvector = particles.mixing_eigenvector

    pbar = tqdm(total = len(particles.list_vtu), desc = "Counting results")
    for i_df in range(len(particles.list_vtu)):

        df = particles.get_df(i_df)

        #print("")

        ID_height = df["ID"][(df["height_color"] > 0) & (df.points[:, 0] < split_height)]
        ID_radius = df["ID"][(df["radius_color"] > 0) & (np.sqrt(np.square(df.points[:, 2]) + np.square(df.points[:, 1])) < split_radius)]
        ID_theta = df["ID"][(df["theta_color"] > 0) & (df.points[:, 2] < 0)]
        
        n_ID_height.append(len(ID_height))
        n_ID_radius.append(len(ID_radius))
        n_ID_theta.append(len(ID_theta))

        if i_df > reference_time_step:
            if n_flows_height < 0 and height_mixing_index[i_df] >= criterion:
                n_flows_height = flow_throughs[i_df]
                time_height = particles.time_list[i_df] - particles.time_list[reference_time_step]

            if n_flows_radius < 0 and radius_mixing_index[i_df] >= criterion:
                n_flows_radius = flow_throughs[i_df]
                time_radius = particles.time_list[i_df] - particles.time_list[reference_time_step]

            if n_flows_theta < 0 and theta_mixing_index[i_df] >= criterion:
                n_flows_theta = flow_throughs[i_df]
                time_theta = particles.time_list[i_df] - particles.time_list[reference_time_step]
            
            if n_flows_doucet < 0 and doucet_mixing_index[i_df] >= criterion:
                n_flows_doucet = flow_throughs[i_df]
                time_doucet = particles.time_list[i_df] - particles.time_list[reference_time_step]

        pbar.update(1)
    
    df_height = pd.DataFrame({"simulation" : i.split('/')[6] + '_' + i.split('/')[5] + '_' + i.split('/')[7].split('_')[-1], "n_flows" : n_flows_height, "time" : time_height}, index=[0])
    df_radius = pd.DataFrame({"simulation" : i.split('/')[6] + '_' + i.split('/')[5] + '_' + i.split('/')[7].split('_')[-1], "n_flows" : n_flows_radius, "time": time_radius}, index=[0])
    df_theta = pd.DataFrame({"simulation" : i.split('/')[6] + '_' + i.split('/')[5] + '_' + i.split('/')[7].split('_')[-1], "n_flows" : n_flows_theta, "time": time_theta}, index=[0])
    df_doucet = pd.DataFrame({"simulation" : i.split('/')[6] + '_' + i.split('/')[5] + '_' + i.split('/')[7].split('_')[-1], "n_flows" : n_flows_doucet, "time": time_doucet}, index=[0])

    if "alginate" in i.split("/") or "Alginate" in i.split("/"):
        folder = "/mnt/d/particles_props/alginate/n_flows.xlsx"
    elif "particle_props_new" in i.split("/"):
        folder = "/mnt/e/particle_props_new/n_flows.xlsx"
    elif "alumina" in i.split("/") or "Alumina" in i.split("/"):
        folder = "/mnt/d/particles_props/alumina/n_flows.xlsx"
    excel = pd.read_excel(folder, sheet_name = None)
    

    excel["height"] = excel["height"].append(df_height, ignore_index = True)
    excel["radius"] = excel["radius"].append(df_radius, ignore_index = True)
    excel["theta"] = excel["theta"].append(df_theta, ignore_index = True)
    excel["doucet"] = excel["doucet"].append(df_doucet, ignore_index = True)

    if drop_duplicates_true:
        excel["height"] = excel["height"].drop_duplicates(subset= "simulation", keep = "last", ignore_index = True)
        excel["radius"] = excel["radius"].drop_duplicates(subset= "simulation", keep = "last", ignore_index = True)
        excel["theta"] = excel["theta"].drop_duplicates(subset= "simulation", keep = "last", ignore_index = True)
        excel["doucet"] = excel["doucet"].drop_duplicates(subset= "simulation", keep = "last", ignore_index = True)

    with pd.ExcelWriter(folder) as writer:
    
        excel["height"].to_excel(writer, sheet_name='height', index = False)
        excel["radius"].to_excel(writer, sheet_name='radius', index = False)
        excel["theta"].to_excel(writer, sheet_name='theta', index = False)
        excel["doucet"].to_excel(writer, sheet_name='doucet', index = False)

    data = np.asarray([particles.time_list[reference_time_step:], flow_throughs[reference_time_step:], n_ID_height[reference_time_step:], n_ID_radius[reference_time_step:], n_ID_theta[reference_time_step:], height_mixing_index[reference_time_step:], radius_mixing_index[reference_time_step:], theta_mixing_index[reference_time_step:], doucet_mixing_index[reference_time_step:], doucet_eigenvector[reference_time_step:][:, 0], doucet_eigenvector[reference_time_step:][:, 1], doucet_eigenvector[reference_time_step:][:, 2], height_mixing_index_std[reference_time_step:], radius_mixing_index_std[reference_time_step:], theta_mixing_index_std[reference_time_step:]])#, doucet_mixing_index[reference_time_step:], doucet_eigenvector[reference_time_step:][:, 0], doucet_eigenvector[reference_time_step:][:, 1], doucet_eigenvector[reference_time_step:][:, 2]

    np.savetxt(f"{i}/data.csv", np.transpose(data), delimiter=',', header = "time, n_flows, n_ID_height, n_ID_radius, n_ID_theta, height_mixing_index, radius_mixing_index, theta_mixing_index,doucet_mixing_index, doucet_eigenvector_r, doucet_eigenvector_theta, doucet_eigenvector_z, height_mixing_index_std, radius_mixing_index_std, theta_mixing_index_std", comments='', fmt='%1.6f')#  doucet_mixing_index, doucet_eigenvector_r, doucet_eigenvector_theta, doucet_eigenvector_z,

    plt.plot(particles.time_list[reference_time_step:], n_ID_height[reference_time_step:], label = "Height")
    plt.plot(particles.time_list[reference_time_step:], n_ID_radius[reference_time_step:], label = "Radius")
    plt.plot(particles.time_list[reference_time_step:], n_ID_theta[reference_time_step:], label = r"$\theta$")
    plt.plot(particles.time_list[reference_time_step:], len(particles.get_df(reference_time_step)["ID"])/4 * doucet_mixing_index[reference_time_step:], label = "Doucet")
    plt.plot(particles.time_list[reference_time_step:], np.repeat(len(particles.get_df(reference_time_step)["ID"])/4, len(particles.time_list[reference_time_step:])), ':', label = "A fourth of particles amount")
    plt.xlabel(r"Time ($s$)")
    plt.ylabel("Number of particles")
    plt.legend()
    plt.savefig(f"{case_path}/{i.split('/')[5] + '_' + i.split('/')[4] + '_' + i.split('/')[6].split('_')[-1]}.png")
    plt.close()

    plt.plot(flow_throughs[reference_time_step:], n_ID_height[reference_time_step:], label = "Height")
    plt.plot(flow_throughs[reference_time_step:], n_ID_radius[reference_time_step:], label = "Radius")
    plt.plot(flow_throughs[reference_time_step:], n_ID_theta[reference_time_step:], label = r"$\theta$")
    plt.plot(flow_throughs[reference_time_step:], len(particles.get_df(reference_time_step)["ID"])/4 * doucet_mixing_index[reference_time_step:], label = "Doucet")
    plt.plot(flow_throughs[reference_time_step:], np.repeat(len(particles.get_df(reference_time_step)["ID"])/4, len(particles.time_list[reference_time_step:])), ':', label = "A fourth of particles amount")
    plt.xlabel(r"Number of flows through")
    plt.ylabel("Number of particles")
    plt.legend()
    plt.savefig(f"{case_path}/{i.split('/')[5] + '_' + i.split('/')[4] + '_' + i.split('/')[6].split('_')[-1]}_flow.png")
    plt.close()

    plt.plot(flow_throughs[reference_time_step:], height_mixing_index[reference_time_step:], label = "Height")
    plt.plot(flow_throughs[reference_time_step:], radius_mixing_index[reference_time_step:], label = "Radius")
    plt.plot(flow_throughs[reference_time_step:], theta_mixing_index[reference_time_step:], label = r"$\theta$")
    plt.plot(flow_throughs[reference_time_step:], doucet_mixing_index[reference_time_step:], label = r"Doucet")
    plt.plot(flow_throughs[reference_time_step:], np.repeat(1, len(particles.time_list[reference_time_step:])), ':', label = "Mixing index = 1")
    plt.xlabel(r"Number of flows through")
    plt.ylabel("Mixing index")
    plt.legend()
    plt.savefig(f"{case_path}/{i.split('/')[5] + '_' + i.split('/')[4] + '_' + i.split('/')[6].split('_')[-1]}_mixing_index.png")
    plt.close()



    #open(subdirectories[i] + "/done.txt", 'a').close()



subdirec = []
subdirectories = listdirs(os.getcwd())

fails = []

if test:
    tempo = time()
    i = subdirectories[test_case_n]
    modify(i)
    print("\n")
    print(f"TIME NEEDED: {(time() - tempo)/60} minutes")
    print("\n")
else:
    for i in subdirectories[first:]:
        try:
            tempo = time()
            modify(i)
            print("\n")
            print(f"FINISHED for: {i}")
            print("\n")
            print(f"TIME NEEDED: {(time() - tempo)/60} minutes")
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