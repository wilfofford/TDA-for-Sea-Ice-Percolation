#  How to use this file:
#
#  make sure the current file is in the same folder as my_functions.py
#
#  set [directory] to be the path to a folder containing 
#  folders named e.g. {Run a1, Run a2,...} where each folder contains files named
#  {Run ai Step j.txt} which store the binary pond images that come from matlab.
#
#  Change [runs] to e.g. ['a0','a1','a2',...]
#  Change [steprange] to e.g. range(1,202) if my timesteps run from 1 to 201
#
#  the program will then create a folder called 'SEDT' containing lots of good things! :)

import numpy as np
import os
import matplotlib.pyplot as plt
import time
from ripser import ripser
from helper_functions import lower_star_distance_matrix, upper_star_distance_matrix_with_disconnected_boundary, SEDTise, persistence_image_bins
from pathlib import Path

directory = '[DIRECTORY PATH]'
runs      = range(1,11)
steprange = range(1,202)


regions = ['A','B','C','D']
ordinal = ['0th','1st','Both']

for r in runs:
    Path(directory+f'/SEDT/0th Homology Sublevel Persistence Images/Run {r}').mkdir(parents=True, exist_ok=True)
    Path(directory+f'/SEDT/1st Homology Sublevel Persistence Images/Run {r}').mkdir(parents=True, exist_ok=True)  
    Path(directory+f'/SEDT/SEDT Numpy/Run {r}').mkdir(parents=True, exist_ok=True)
    Path(directory+f'/SEDT/Persistence Diagrams/Run {r}').mkdir(parents=True, exist_ok=True)
    Path(directory+f'/SEDT/Visualization/Run {r}').mkdir(parents=True, exist_ok=True)
    
    for o in ordinal:
        Path(directory+f'/SEDT/For Bokeh/{o} Hom/Run {r}').mkdir(parents=True, exist_ok=True)
    
    for rg in regions:        
        Path(directory+f'/SEDT/Regions/{rg}/Run {r}').mkdir(parents=True, exist_ok=True)
        
    Path(directory+f'/SEDT/Superlevel Region A (With Boundary)/Run {r}').mkdir(parents=True, exist_ok=True)


ponds={}
for r in runs:
    for s in steprange:
        ponds[(r,s)]=np.loadtxt(directory+f'/Run {r}/Run {r} Step {s}.txt', delimiter=',')




SEDT={}
sublevel={}
superlevel={}




for r,s in ponds:
    print(f'working on: run {r} step {s} (calculating homology)')
    start=time.perf_counter()
          
    SEDT[(r,s)]=SEDTise(ponds[r,s])
    
    os.chdir(directory+f'/SEDT/SEDT Numpy/Run {r}')
    np.save(f'Run {r} Step {s} SEDT.npy', SEDT[r,s])
          
    sublevel[(r,s)]=ripser(lower_star_distance_matrix(SEDT[r,s]), distance_matrix=True)['dgms']
    os.chdir(directory+f'/SEDT/Persistence Diagrams/Run {r}')
    np.save(f'Run {r} Step {s} Persistence Diagram.npy', sublevel[r,s])
    


    superlevel[(r,s)]=ripser(upper_star_distance_matrix_with_disconnected_boundary(SEDT[r,s]), distance_matrix=True)['dgms']

    A=superlevel[r,s][0][np.logical_and(superlevel[r,s][0][:,0]<0,superlevel[r,s][0][:,1]>0)]
    os.chdir(directory+f'/SEDT/Superlevel Region A (With Boundary)/Run {r}')
    np.save(f'Run {r} Step {s} Superlevel Region A (With Boundary).npy',A)
          
          
    end=time.perf_counter()
    print(f'runtime: {end-start}s')
    



#                       Regions:
#  We break the persistence diagrams for the SEDT (sublevel) into 4 regions:
#  
#  0th Homology:                1st Homology:
#
#       |                            |
#   'A' |                        'C' | 'D'
#   ----+----                    ----+----
#   'B' |                            |
#       |                            |
#
#  These regions reflect various aspect of the melt pond geometry.
#  This also explains why we use a 'binning' method for constructing persistence
#  images; we don't want to use smooth gaussian kernels because they will interfere
#  between different regions


for key in sublevel:
    A=sublevel[key][0][np.logical_and(sublevel[key][0][:,0]<0,sublevel[key][0][:,1]>0)]
    os.chdir(directory + f'/SEDT/Regions/A/Run {key[0]}')
    np.save(f'Run {key[0]} Step {key[1]} Region A',A)
    
    B=sublevel[key][0][np.logical_and(sublevel[key][0][:,0]<0,sublevel[key][0][:,1]<0)]
    os.chdir(directory+f'/SEDT/Regions/B/Run {key[0]}')
    np.save(f'Run {key[0]} Step {key[1]} Region B',B)
    
    C=sublevel[key][1][np.logical_and(sublevel[key][1][:,0]<0,sublevel[key][1][:,1]>0)]
    os.chdir(directory+f'/SEDT/Regions/C/Run {key[0]}')
    np.save(f'Run {key[0]} Step {key[1]} Region C',C)
    
    D=sublevel[key][1][np.logical_and(sublevel[key][1][:,0]>0,sublevel[key][1][:,1]>0)]
    os.chdir(directory+f'/SEDT/Regions/D/Run {key[0]}')
    np.save(f'Run {key[0]} Step {key[1]} Region D',D)
    
#this function lives in my_functions.py , where I've documented what each argument does

pimages = persistence_image_bins(sublevel, save=True, visualize=True, SEDTs=SEDT, Ponds=ponds, 
        directoryfunction = lambda key : (directory+f'/SEDT/0th Homology Sublevel Persistence Images/Run {key[0]}',
                                          directory+f'/SEDT/1st Homology Sublevel Persistence Images/Run {key[0]}',
                                          directory+f'/SEDT/Visualization/Run {key[0]}',
                                          directory+f'/SEDT/For Bokeh/0th Hom/Run {key[0]}',
                                          directory+f'/SEDT/For Bokeh/1st Hom/Run {key[0]}',
                                          directory+f'/SEDT/For Bokeh/Both Hom/Run {key[0]}'),
        filenamefunction  = lambda key : f'Run {key[0]} Step {key[1]}',
        plottitlefunction = lambda key : f'Run {key[0]} Step {key[1]}')
          

