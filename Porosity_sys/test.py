list1 = [0.0, 4.0, 8.0, 6.5, 10.5, 2.5, 5.0, 9.0, 1.0, 7.5, 3.5, 2.0, 10.0, 6.0, 0.5, 4.5, 8.5, 7.0, 11.0, 3.0, 5.5, 9.5, 1.5]
list2 = ['CFD_0.vtk', 'CFD_34000.vtk', 'CFD_10000.vtk', 'CFD_6000.vtk', 'CFD_30000.vtk', 'CFD_12000.vtk', 'CFD_38000.vtk', 'CFD_14000.vtk', 'CFD_16000.vtk', 'CFD_42000.vtk', 'CFD_18000.vtk', 'CFD_44000.vtk', 'CFD_32000.vtk', 'CFD_2000.vtk', 'CFD_8000.vtk', 'CFD_20000.vtk', 'CFD_22000.vtk', 'CFD_36000.vtk', 'CFD_24000.vtk', 'CFD_40000.vtk', 'CFD_28000.vtk', 'CFD_26000.vtk', 'CFD_4000.vtk']
list1, list2 = zip(*sorted(zip(list1, list2)))

print(list1)
print(list2)