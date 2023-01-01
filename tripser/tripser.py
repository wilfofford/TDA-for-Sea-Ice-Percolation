import numpy as np
import math
import os
import subprocess

"""
tripser - ripser with TRacking. computes 0th homology barcodes for lower star filtration
and records the pixel coordinates of birth/death events
input: 2D array of floats
output: dictionary with 'dgm': standard 0th homology barcode
                       'locations': Nx4 array with rows [birth y_coord, birth x_coord, death y_coord, death]
                       
make sure compiled c++ program tripser.exe is in same directory
"""
                              
def tripser(img):
    m,n=img.shape
    flat=img.flatten().tolist()
    dgm=[]
    locations=[]

    with open("storing_img.txt","w") as f:
        f.write(str(m)+"\n"+str(n))
        for item in flat:
            f.write("\n"+str(item))

    process = subprocess.Popen(["./tripser"])
    process.wait()
    

    with open("storing_barcode.txt","r") as g:
        for line in g:
            chunks=line.split(' ')
            dgm.append([float(chunks[i]) for i in range(2)])
            locations.append([int(chunks[i]) for i in range(2,6)])
    
    os.remove("storing_img.txt")
    os.remove("storing_barcode.txt")
            
    return{'dgm':np.array(dgm),'locations':np.array(locations)}

def compute_distances(barcode):
    N=barcode['dgm'].shape[0]
    distDGM=np.zeros([N,3])
    distDGM[:,:2]=barcode['dgm']
    for i in range(N):
        distDGM[i,2]=math.sqrt((barcode['locations'][i][0]-barcode['locations'][i][2])**2+(barcode['locations'][i][1]-barcode['locations'][i][3])**2)

    return distDGM