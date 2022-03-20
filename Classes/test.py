from Lethe_pyvista_tools import *

first_case = Lethe_pyvista_tools()
first_case.lethe_import_to_pyvista('result_particles.pvd')

print(first_case.list_vtu)