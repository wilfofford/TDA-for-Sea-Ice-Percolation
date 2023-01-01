import numpy as np
import math
import copy
from scipy import sparse
from functools import reduce
from scipy.interpolate import griddata
from scipy.ndimage.morphology import distance_transform_edt as edt
import warnings
from persim import plot_diagrams
import os
import matplotlib.pyplot as plt
import time
import matplotlib.colors
pond_cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["#a1bcc1","#153440"])



#helper functions used in other files




#takes in a binary (true/false) image and outputs the signed euclidean distance transform. 'true' values are
#sent to negative values, 'false' are sent to positive.
def SEDTise(img):
    pass1=edt(img)
    pass2=edt(1-img)
    
    return pass2-pass1

# def SEDTise(img):                                             

#     SEDT=copy.deepcopy(img).astype(float)

        
#     for j in range(img.shape[0]):
#         for k in range(img.shape[1]):
#             if not np.isnan(img[j,k]):
#                 if img[j,k]:
#                     searchradius = 1
                
#                     while img[max(0,j-searchradius):min(img.shape[0]-1, j+searchradius)+1,
#                                                 max(0,k-searchradius):min(img.shape[1]-1,k+searchradius)+1].all():
#                         searchradius+=1
                        
#                     subsquare=img[max(0,j-searchradius):min(img.shape[0]-1, j+searchradius)+1,
#                                                 max(0,k-searchradius):min(img.shape[1]-1,k+searchradius)+1]
#                     dists=[]
                    
                    
#                     for r in range(subsquare.shape[0]):
#                         for l in range(subsquare.shape[1]):
#                             if not subsquare[r,l]:
#                                 dists.append(math.sqrt((r+max(0,j-searchradius)-j)**2+(l+max(0,k-searchradius)-k)**2))
                                
                    
#                     SEDT[j,k]=-min(dists)
                    
#                 if not img[j,k]:
#                     searchradius = 1
                
#                     while not (img[max(0,j-searchradius):min(img.shape[0]-1, j+searchradius)+1,
#                                                 max(0,k-searchradius):min(img.shape[1]-1,k+searchradius)+1]).any():
#                         searchradius+=1
                        
#                     subsquare=img[max(0,j-searchradius):min(img.shape[0]-1, j+searchradius)+1,
#                                                 max(0,k-searchradius):min(img.shape[1]-1,k+searchradius)+1]
#                     dists=[]
                    
#                     for r in range(subsquare.shape[0]):
#                         for l in range(subsquare.shape[1]):
#                             if subsquare[r,l]:
#                                 dists.append(math.sqrt((r+max(0,j-searchradius)-j)**2+(l+max(0,k-searchradius)-k)**2))
                                
#                     SEDT[j,k]=min(dists)
#     return SEDT


def lower_star_distance_matrix(img):
    """
    Construct a lower star filtration on an image
    Parameters
    ----------
    img: ndarray (M, N)
        An array of single channel image data
    Returns
    -------
    I: ndarray (K, 2)
        A 0-dimensional persistence diagram corresponding to the sublevelset filtration
    """
    m, n = img.shape
    idxs = np.arange(m * n).reshape((m, n))
    I = idxs.flatten()
    J = idxs.flatten()
    V = img.flatten()
    # Connect 8 spatial neighbors
    tidxs = np.ones((m + 2, n + 2), dtype=np.int64) * np.nan
    tidxs[1:-1, 1:-1] = idxs
    tD = np.ones_like(tidxs) * np.nan
    tD[1:-1, 1:-1] = img
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            if di == 0 and dj == 0:
                continue
            thisJ = np.roll(np.roll(tidxs, di, axis=0), dj, axis=1)
            thisD = np.roll(np.roll(tD, di, axis=0), dj, axis=1)
            thisD = np.maximum(thisD, tD)
            # Deal with boundaries
            boundary = ~np.isnan(thisD)
            thisI = tidxs[boundary]
            thisJ = thisJ[boundary]
            thisD = thisD[boundary]
            I = np.concatenate((I, thisI.flatten()))
            J = np.concatenate((J, thisJ.flatten()))
            V = np.concatenate((V, thisD.flatten()))
    return sparse.coo_matrix((V, (I, J)), shape=(idxs.size, idxs.size))

def upper_star_distance_matrix_with_disconnected_boundary(img0):
    #This function returns an upper star distance matrix, with the following alterations:
    #Pixels are now 4-connected, not 8-connected
    #A boundary is added, but each boundary pixel is not connected to other boundary pixels.
    img=-img0
    m, n = img.shape
    idxs = np.arange(m * n).reshape((m, n))
    I = idxs.flatten()
    J = idxs.flatten()
    V = img.flatten()
    # Connect 4 spatial neighbors
    tidxs = np.ones((m + 2, n + 2), dtype=np.int64) * np.nan
    tidxs[1:-1, 1:-1] = idxs
    tD = np.ones_like(tidxs) * np.nan
    tD[1:-1, 1:-1] = img
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            #don't connect diagonally 
            if np.abs(di)+np.abs(dj)==2:
                continue
            if di == 0 and dj == 0:
                continue
            thisJ = np.roll(np.roll(tidxs, di, axis=0), dj, axis=1)
            thisD = np.roll(np.roll(tD, di, axis=0), dj, axis=1)
            thisD = np.maximum(thisD, tD)
            # Deal with boundaries
            boundary = ~np.isnan(thisD)
            thisI = tidxs[boundary]
            thisJ = thisJ[boundary]
            thisD = thisD[boundary]
            I = np.concatenate((I, thisI.flatten()))
            
            J = np.concatenate((J, thisJ.flatten()))
            
            V = np.concatenate((V, thisD.flatten()))
    
    
    K=m*n
    #add disconnected boundary with filtration value -inf
    outer_elements=np.arange(m)
    for i in range(n-2):
       
        outer_elements = np.concatenate((outer_elements,np.array([m*(i+1),m*(i+2)-1])))
    outer_elements=np.concatenate((outer_elements,np.arange(m*(n-1),m*n)))
                                            
    L=outer_elements.size
    boundary_coordinates=np.arange(K,K+L)
    
    boundary_birthtimes=np.zeros(L)
    for i in range(L):
        
        #ice components on the boundary are born from the beginning.
        #when a component that did not start off touching the boundary connects to the boundary, it dies.
        if V[outer_elements[i]]<0:
            boundary_birthtimes[i]=-np.inf#V[outer_elements[i]]
        else:
            boundary_birthtimes[i]=-1000000  #ideally the first -np.inf would be "-2 infinity" and this -1million would be "-1 infinity"
    
    
    
    I=np.concatenate((I, boundary_coordinates, boundary_coordinates, outer_elements))
    J=np.concatenate((J, boundary_coordinates, outer_elements, boundary_coordinates))
    V=np.concatenate((V, boundary_birthtimes, V[outer_elements], V[outer_elements]))
    
    return sparse.coo_matrix((V, (I, J)), shape=(K+L, K+L))

# def upper_star_distance_matrix(img0):
#     """
#     Construct a lower star filtration on an image
#     Parameters
#     ----------
#     img: ndarray (M, N)
#         An array of single channel image data
#     Returns
#     -------
#     I: ndarray (K, 2)
#         A 0-dimensional persistence diagram corresponding to the sublevelset filtration
#     """
#     #I've just altered one line to stop it recording diagonal connections  
#     img=-img0
#     m, n = img.shape
#     idxs = np.arange(m * n).reshape((m, n))
#     I = idxs.flatten()
#     J = idxs.flatten()
#     V = img.flatten()
#     # Connect 8 spatial neighbors
#     tidxs = np.ones((m + 2, n + 2), dtype=np.int64) * np.nan
#     tidxs[1:-1, 1:-1] = idxs
#     tD = np.ones_like(tidxs) * np.nan
#     tD[1:-1, 1:-1] = img
#     for di in [-1, 0, 1]:
#         for dj in [-1, 0, 1]:
#             #don't connect diagonally 
#             if np.abs(di)+np.abs(dj)==2:
#                 continue
#             if di == 0 and dj == 0:
#                 continue
#             thisJ = np.roll(np.roll(tidxs, di, axis=0), dj, axis=1)
#             thisD = np.roll(np.roll(tD, di, axis=0), dj, axis=1)
#             thisD = np.maximum(thisD, tD)
#             # Deal with boundaries
#             boundary = ~np.isnan(thisD)
#             thisI = tidxs[boundary]
#             thisJ = thisJ[boundary]
#             thisD = thisD[boundary]
#             I = np.concatenate((I, thisI.flatten()))
#             J = np.concatenate((J, thisJ.flatten()))
#             V = np.concatenate((V, thisD.flatten()))
#     return sparse.coo_matrix((V, (I, J)), shape=(idxs.size, idxs.size))


def interp(raw):

    x = np.arange(0, raw.shape[1])
    y = np.arange(0, raw.shape[0])
    #mask invalid values
    array = np.ma.masked_invalid(raw)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]

    interpolated = griddata((x1, y1), newarr.ravel(),
                            (xx, yy),
                                method='linear')
    
    #other interpolation methods: 'nearest', 'cubic'
    
    return interpolated

def kNNimg(img, k, pixelwidth=1):
    #takes in binary image of points (nan is counted as false), returns codensity (sum of distances to kth nearest neigbours)
    
    
    codensityOutput=np.zeros((img.shape[0],img.shape[1],2))+np.nan
    codensityOutput[:,:,0]=img[:,:]
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            
            
            if img[i,j] and not np.isnan(img[i,j]):
                searchradius=0
                count=1
                while count<k+1:
                    searchradius+=1
                    toplefty=max(0,i-searchradius)
                    topleftx=max(0,j-searchradius)
                    bottomrighty=min(img.shape[0]-1,i+searchradius)
                    bottomrightx=min(img.shape[1]-1,j+searchradius)
                    
                    
                    
                    count=(np.logical_and(img[toplefty:bottomrighty+1,topleftx:bottomrightx+1], 
                                          ~np.isnan(img[toplefty:bottomrighty+1,topleftx:bottomrightx+1]))).sum()
                
                
                toplefty=max(0,i-searchradius)
                topleftx=max(0,j-searchradius)
                bottomrighty=min(img.shape[0]-1,i+searchradius)
                bottomrightx=min(img.shape[1]-1,j+searchradius)
                xrange=bottomrightx-topleftx+1
                yrange=bottomrighty-toplefty+1
                
                
                distances=np.zeros(xrange*yrange)+np.nan
                for m in range(toplefty,bottomrighty+1):
                    for n in range(topleftx,bottomrightx+1):
                        if img[m,n] and not np.isnan(img[m,n]):
                            m0=m-toplefty
                            n0=n-topleftx
                            distances[n0*yrange+m0]=math.sqrt(((i-m)*pixelwidth)**2+((j-n)*pixelwidth)**2)
                distancesNotNan=distances[~np.isnan(distances)]
                
                codensityOutput[i,j,1]=np.sum(np.sort(distancesNotNan)[:k+1])
    return codensityOutput






#  takes in a dictionary of persistence diagrams, and outputs a dictionary containing {key:(0th homology persistence 
#  image, 1st hom Persistence Image)}. persistence images are computed by 'bins' of size 1.
#  if save=True, then the function will save these.
#  if visualize=True, then it will also save plots of the diagrams, SEDT and images

#  these will be saved to the directory specified by directoryfunction, e.g.
#  lambda key : ('0th Hom directory', '1st Hom Directory', *optional plots directory*, *optional thumbnail directory for 0th hom*,
#  *optional thumbnail directory for 1st hom*, *optional thumbnail directory for both hom*)
#  where key ranges through the keys of dgms

#  similarly, the file name is determined by filename function e.g.
#  lambda key : f'run {key[0]} step {key[1]}'

#  note: if visualize=True, save must be true and you must enter Ponds and SEDTs

def persistence_image_bins(dgms, save=False, visualize=False, directoryfunction=None, filenamefunction=None, Ponds=None, 
SEDTs=None, plottitlefunction=None):

    if visualize and not save:
        warnings.warn('You have set visualize to true, but save is set to false, so visualize will do nothing. Please set save to true and try again.')

    if visualize and (Ponds == None or SEDTs == None or directoryfunction == None or filenamefunction == None or plottitlefunction == None):
        warnings.warn('You have set visualize to true but not entered all SEDT images, Pond images and directory functions')

    
    min0=math.floor(min([np.nanmin(dgms[key][0]) for key in dgms]))
    max0=math.ceil(max([np.nanmax(dgms[key][0][dgms[key][0]!=np.inf]) for key in dgms]))
    min1=math.floor(min([np.nanmin(dgms[key][1]) for key in dgms]))
    max1=math.ceil(max([np.nanmax(dgms[key][1][dgms[key][1]!=np.inf]) for key in dgms]))

    diff0=max0-min0
    diff1=max1-min1
    

    pimages0={key:np.zeros((diff0+1,diff0+1)) for key in dgms}
    pimages1={key:np.zeros((diff1+1,diff1+1)) for key in dgms}



    for key in dgms:
        for [b,d] in dgms[key][0]:
                if not np.isnan(b) and not np.isnan(d) and d!=np.inf:
                    pimages0[key][max0-int(np.floor(d)), int(np.floor(b))-min0]+=1


        for [b,d] in dgms[key][1]:
                if not np.isnan(b) and not np.isnan(d) and d!=np.inf:
                    pimages1[key][max1-int(np.floor(d)), int(np.floor(b))-min1]+=1

    output={key:(pimages0[key],pimages1[key]) for key in dgms}
    
    
    
    if save:
        for key in dgms:
            os.chdir(directoryfunction(key)[0])
            np.save(filenamefunction(key) + ' 0th Homology Persistence Image.npy', pimages0[key])
            
            os.chdir(directoryfunction(key)[1])
            np.save(filenamefunction(key) + ' 1st Homology Persistence Image.npy', pimages1[key])
            
        if visualize:
            
            for key in dgms:
                start=time.perf_counter()
                print(f'working on: Run {key[0]} Step {key[1]} (plotting images)')

                os.chdir(directoryfunction(key)[2])
            
                plt.figure(figsize=(30,20))

                plt.subplot(231)
                plt.imshow(Ponds[key], cmap=pond_cmap)
                plt.title (plottitlefunction(key))

                plt.subplot(234)
                plt.imshow(SEDTs[key])

                plt.subplot(232)
                plt.imshow(np.log(pimages0[key]+1))
                plt.plot([0,diff0],[diff0,0], lw=0.3, c='white')
                plt.plot([-min0,-min0], [diff0,0], lw=0.3, c='white')
                plt.plot([0,diff0],[max0,max0], lw=0.3, c='white')

                plt.title('0th Homology')

                plt.subplot(235)
                plot_diagrams(dgms[key][0], labels='$H_0$', xy_range=[min0,max0]*2, size=20)

                plt.subplot(233)
                plt.imshow(np.log(pimages1[key]+1))
                plt.plot([0,diff1],[diff1,0], lw=0.3, c='white')
                plt.plot([-min1,-min1], [diff1,0], lw=0.3, c='white')
                plt.plot([0,diff1],[max1,max1], lw=0.3, c='white')

                plt.title('1st Homology')

                plt.subplot(236)
                plot_diagrams(dgms[key][1], labels='$H_1$', xy_range=[min1,max1]*2, size=20)

                os.chdir(directoryfunction(key)[2])
                plt.savefig(filenamefunction(key) + ' Plot.png')

                plt.close('all')
                end=time.perf_counter()
                print(f'runtime: {end-start}s')

            #  We also save some smaller thumbnails that we can use for visualizing PCA

            
            for key in dgms:
                start = time.perf_counter()
                print(f'working on: run {key[0]} step {key[1]} (plotting thumbnails)')

                os.chdir(directoryfunction(key)[3])
            
                plt.figure(figsize=(20,10))
                plt.subplot(121)
                plt.imshow(Ponds[key], cmap=pond_cmap)
                plt.title(f'Run {key[0]} Step {key[1]}')

                plt.subplot(122)
                plt.imshow(np.log(pimages0[key]+1))
                plt.plot([0,diff0],[diff0,0], lw=0.3, c='white')
                plt.plot([-min0,-min0], [diff0,0], lw=0.3, c='white')
                plt.plot([0,diff0],[max0,max0], lw=0.3, c='white')
                plt.title('0th Homology')

                plt.savefig(filenamefunction(key)+' 0th Homology Thumbnail.png')
                
                plt.close('all')


                os.chdir(directoryfunction(key)[4])

                plt.figure(figsize=(20,10))

                plt.subplot(121)
                plt.imshow(Ponds[key], cmap=pond_cmap)
                plt.title(f'Run {key[0]} Step {key[1]}')

                plt.subplot(122)
                plt.imshow(np.log(pimages1[key]+1))
                plt.plot([0,diff1],[diff1,0], lw=0.3, c='white')
                plt.plot([-min1,-min1], [diff1,0], lw=0.3, c='white')
                plt.plot([0,diff1],[max1,max1], lw=0.3, c='white')
                plt.title(f'1st Homology')

                plt.savefig(filenamefunction(key)+' 1st Homology Thumbnail.png')

                plt.close('all')


                os.chdir(directoryfunction(key)[5])

                plt.figure(figsize=(30,10))
                plt.subplot(131)
                plt.imshow(Ponds[key], cmap=pond_cmap)
                plt.title(f'Run {key[0]} Step {key[1]}')

                plt.subplot(132)
                plt.imshow(np.log(pimages0[key]+1))
                plt.plot([0,diff0],[diff0,0], lw=0.3, c='white')
                plt.plot([-min0,-min0], [diff0,0], lw=0.3, c='white')
                plt.plot([0,diff0],[max0,max0], lw=0.3, c='white')
                plt.title('0th Homology')

                plt.subplot(133)
                plt.imshow(np.log(pimages1[key]+1))
                plt.plot([0,diff1],[diff1,0], lw=0.3, c='white')
                plt.plot([-min1,-min1], [diff1,0], lw=0.3, c='white')
                plt.plot([0,diff1],[max1,max1], lw=0.3, c='white')
                plt.title(f'1st Homology')

                plt.savefig(filenamefunction(key)+' Both Homology Thumbnail.png')

                plt.close('all')

                end=time.perf_counter()
                print(f'runtime: {end-start}s')
                
    return output


def emptymax(arr):
    if arr.size:
        return np.nanmax(arr)
    else:
        return np.nan

def emptymin(arr):
    if arr.size:
        return np.nanmin(arr)
    else:
        return np.nan



#truncates values if they are below a certain proportion of the data range
def truncated_mean(arr, truncate_proportion):
    M=emptymax(arr)
    m=emptymin(arr)
    threshold=m+truncate_proportion*(M-m)

    list1=arr.tolist()
    list_truncated=[item for item in list1 if not np.isnan(item) and not np.isinf(item) and item>=threshold]

    return np.mean(list_truncated)



#returns boundary of binary image
def binary_boundary(img0, upscale=4):
    I0=img0.shape[0]
    J0=img0.shape[1]
    img=np.zeros((I0*upscale,J0*upscale))
    for i in range(I0):
        for j in range(J0):
            img[upscale*i:upscale*i+upscale,upscale*j:upscale*j+upscale]=img0[i,j]
    
    I=img.shape[0]
    J=img.shape[1]
    boundary=np.zeros_like(img)
    for i in range(I):
        for j in range(J):
            min_v=max(0,i-1)
            min_h=max(0,j-1)
            max_v=min(i+1,I-1)
            max_h=min(j+1,J-1)
            
            if not img[min_v:max_v+1,min_h:max_h+1].all() and (img[min_v:max_v+1,min_h:max_h+1]).any():
                boundary[i,j]=1
    return boundary
    
