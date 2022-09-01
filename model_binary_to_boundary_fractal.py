#  How to use this file:
#
#  set [directory] to be the path to a folder containing 
#  folders named e.g. {Run 1, Run 2,...} where each folder contains files named
#  {Run i Step j.txt} which store the binary pond images that come from matlab.
#
#  Change [runs] to e.g. [0,1,2,...]
#  Change [steprange] to e.g. range(1,202) if my timesteps run from 1 to 201
#
#  the program will then create a folder called 'Fractal of Boundary' containing lots of good things! :)





import numpy as np
import os
import matplotlib.pyplot as plt
import time
from ripser import ripser
from pathlib import Path
import random
from scipy.stats import linregress
import pandas as pd
from my_functions import binary_boundary



def ran(first,last,num_steps):
   increment0=(last-first)/num_steps
   list0 = [first]
   val=first
   for i in range(num_steps):
       val+=increment0
       list0.append(val)
   return(list0)


directory = 'DIRECTORY_PATH'
runs      = range(1,11)
steprange = range(1,202)
do_1st_hom=False

#maximum points to subsample from image.
max_points=10000
alphas=[0.25,0.5,0.75,1,1.25,1.5]

pond_image={}
pond_boundary={}
pointcloud ={}
E0_graph={}
log_E0_graph={}

maxdim=0

if do_1st_hom:
    E1_graph={}
    log_E1_graph={}
    maxdim=1


#create folders 
for r in runs:
    for alpha in alphas:
        Path(directory+f'/Fractal of Boundary/log E0 graph alpha {alpha}/Run {r}').mkdir(parents=True, exist_ok=True)
        Path(directory+f'/Fractal of Boundary/Regression Plots E0/Run {r}').mkdir(parents=True, exist_ok=True)
        
        if do_1st_hom:
            Path(directory+f'/Fractal of Boundary/log E1 graph alpha {alpha}/Run {r}').mkdir(parents=True, exist_ok=True)
            Path(directory+f'/Fractal of Boundary/Regression Plots E1/Run {r}').mkdir(parents=True, exist_ok=True)

#load files
for r in runs:
    for s in steprange:
        pond_image[(r,s)]=np.loadtxt(directory+f'/Run {r}/Run {r} Step {s}.txt', delimiter=',')
        pond_boundary[(r,s)]=binary_boundary(pond_image[r,s])
        
        
#convert image to point cloud
for key in pond_boundary:
    I = pond_boundary[key].shape[0]
    J = pond_boundary[key].shape[1]
    flatlist=[]
    for i in range(I):
        for j in range(J):
            if pond_boundary[key][i,j]:
                flatlist.append([j,i])
    pointcloud[key]=np.array(flatlist)

for r in runs:
    for s in steprange:
        print(f'Working on: Run {r} Step {s}')
        start0=time.perf_counter()
        
        for a in alphas:
            log_E0_graph[(a,r,s)]=np.zeros((25,2))
            E0_graph[(a,r,s)]=np.zeros((25,2))
            
            if do_1st_hom:
                log_E1_graph[(a,r,s)]=np.zeros((25,2))
                E1_graph[(a,r,s)]=np.zeros((25,2))

        N=pointcloud[r,s].shape[0]
        big=min(max_points, max(N,500))
        
        #choose logarithmically spaced sample sizes
        values =[int(np.floor(np.exp(j))) for j in ran(np.log(big/2),np.log(big),24)]

        for m in range(25):
            print(f'Number of points: {values[m]}')
            start1=time.perf_counter()
            sample=np.zeros((values[m],2))
            for k in range(values[m]):
                
               #choose random sample, with continuous uniform noise 
               sample[k]=pointcloud[r,s][random.randrange(N)]+[random.random()-0.5,random.random()-0.5]
            diagram=ripser(sample,maxdim=maxdim, distance_matrix=False)['dgms']

            
            #create list of persistence values
            diagram0=diagram[0][:-1]
            Ivalues0=np.zeros(diagram0.shape[0])
            for t in range(diagram0.shape[0]):
                Ivalues0[t]=((diagram0[t,1]-diagram0[t,0]))
            E0={}
            
            if do_1st_hom:    
                diagram1=diagram[1]
                Ivalues1=np.zeros(diagram1.shape[0])
                for t in range(diagram1.shape[0]):
                    Ivalues1[t]=((diagram1[t,1]-diagram1[t,0]))
                E1={}

            #calculate Ei statistic i.e. sum of persistence^alpha
            for a in alphas:
                E0[a]=np.sum(Ivalues0**a)
                E0_graph[a,r,s][m]=[values[m],E0[a]]
                
                if do_1st_hom:
                    E1[a]=np.sum(Ivalues1**a)
                    E1_graph[a,r,s][m]=[values[m],E1[a]]
            
            end1=time.perf_counter()
            print(f'Runtime: {end1-start1}s')
        
        #save log/log graph
        for a in alphas:
            log_E0_graph[a,r,s]=np.log(E0_graph[a,r,s])
            os.chdir(directory+f'/Fractal of Boundary/log E0 graph alpha {a}/Run {r}')
            np.save(f'Run {r} Step {s} log E0 graph alpha {a}.npy',log_E0_graph[a,r,s])

            if do_1st_hom:
                log_E1_graph[a,r,s]=np.log(E1_graph[a,r,s])
                os.chdir(directory+f'/Fractal of Boundary/log E1 graph alpha {a}/Run {r}')
                np.save(f'Run {r} Step {s} log E1 graph alpha {a}.npy',log_E1_graph[a,r,s])


        end0=time.perf_counter()
        print(f'Total time: {end0-start0}s')

colours={0.25:'coral', 0.5:'gold', 0.75:'yellowgreen', 1:'darkturquoise', 1.25:'mediumpurple', 1.5:'orchid'}

xx0={}
yy0={}
slopes0={}
intercepts0={}
dim0={}

if do_1st_hom:
    slopes1={}
    intercepts1={}
    xx1={}
    yy1={}
    dim1={}

#plot regression
for r in runs:
    for s in steprange:
        print(f'Plotting Run {r} Step {s} Regression')
        for a in alphas:
            xx0[(a,r,s)]=log_E0_graph[a,r,s][:,0]
            yy0[(a,r,s)]=log_E0_graph[a,r,s][:,1]
            slopes0[(a,r,s)]=linregress(xx0[a,r,s],yy0[a,r,s]).slope
            intercepts0[(a,r,s)]=linregress(xx0[a,r,s],yy0[a,r,s]).intercept
            dim0[(a,r,s)]=a/(1-slopes0[a,r,s])
            
            
            if do_1st_hom:
                xx1[(a,r,s)]=log_E1_graph[a,r,s][:,0]
                yy1[(a,r,s)]=log_E1_graph[a,r,s][:,1]
                slopes1[(a,r,s)]=linregress(xx1[a,r,s],yy1[a,r,s]).slope
                intercepts1[(a,r,s)]=linregress(xx1[a,r,s],yy1[a,r,s]).intercept            
                dim1[(a,r,s)]=a/(1-slopes1[a,r,s])

        fig, ax=plt.subplots()
        plt.title(f'Run {r} Step {s} 0th Homology Regression')

        for a in alphas:
            ax.scatter(xx0[a,r,s],yy0[a,r,s],
            label=f'α = {a}, dimension = {round(dim0[a,r,s],4)}',
            c=colours[a])

            x1=np.min(xx0[a,r,s])
            x2=np.max(xx0[a,r,s])
            y1=slopes0[a,r,s]*x1+intercepts0[a,r,s]
            y2=slopes0[a,r,s]*x2+intercepts0[a,r,s]

            ax.plot([x1,x2],[y1,y2],c=colours[a])
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*0.8,box.height])

        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))

        os.chdir(directory+f'/Fractal of Boundary/Regression Plots E0/Run {r}')
        plt.savefig(f'Run {r} Step {s} 0th Homology Regression.png')

        plt.close('all')



        if do_1st_hom:
            fig, ax=plt.subplots()
            plt.title(f'Run {r} Step {s} 1st Homology Regression')

            for a in alphas:
                ax.scatter(xx1[a,r,s],yy1[a,r,s],
                label=f'α = {a}, dimension = {round(dim1[a,r,s],4)}',
                c=colours[a])

                x1=np.min(xx1[a,r,s])
                x2=np.max(xx1[a,r,s])
                y1=slopes1[a,r,s]*x1+intercepts1[a,r,s]
                y2=slopes1[a,r,s]*x2+intercepts1[a,r,s]

                ax.plot([x1,x2],[y1,y2],c=colours[a])
            box = ax.get_position()
            ax.set_position([box.x0,box.y0,box.width*0.8,box.height])

            ax.legend(loc='center left', bbox_to_anchor=(1,0.5))

            os.chdir(directory+f'/Fractal of Boundary/Regression Plots E1/Run {r}')
            plt.savefig(f'Run {r} Step {s} 1st Homology Regression.png')

            plt.close('all')

os.chdir(directory+'/Fractal')

steps=list(steprange)

fractal_dim=pd.DataFrame(steps, columns=['Step'])

#save fractal dimension dataframe
for r in runs:
    for a in alphas:
        fractal_dim[f'Run {r} alpha {a} PH0 dim']=[dim0[a,r,s] for s in steprange]
        
        if do_1st_hom:
            fractal_dim[f'Run {r} alpha {a} PH1 dim']=[dim1[a,r,s] for s in steprange]
fractal_dim.to_pickle(directory+'/Fractal of Boundary/Boundary Fractal Dimensions.pkl')
