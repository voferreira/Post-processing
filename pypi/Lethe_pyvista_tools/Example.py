from Lethe_pyvista_tools import *

#This script prints the content of your prm file as dictionary
#To run Lethe_pyvista_tools you need to specify the path to your
#case and the name of the .prm file
example = Lethe_pyvista_tools('PATH TO YOUR CASE', 'NAME_OF_YOUR_PARAMETERS_FILE.prm')

print('This is the dictionary of your .prm file:')
print(example.prm_dict)
print('To print out any value in the dictionary, ask for it using ["parameter_name"] right after .prm_dict variable')

#To read the data to pyvista dataframe, use the following with
#the .pvd file as argument

example.read_lethe_to_pyvista('NAME_OF_YOUR_PVD_FILE.pvd')