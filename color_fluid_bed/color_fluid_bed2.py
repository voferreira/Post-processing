

import matplotlib.pyplot as plt

from lethe_pyvista_tools import *

case_path = "/mnt/d/particles_props/restitution/high_velocity/liquid_fluidized_bed_0.01"

particles = lethe_pyvista_tools(case_path = case_path, prm_file_name = "liquid_fluidized_bed.prm")

particles.read_lethe_to_pyvista(pvd_name = "mod_result__particles.pvd")

particles.mixing_index_nearest_neighbors()

print(particles.mixing_index)