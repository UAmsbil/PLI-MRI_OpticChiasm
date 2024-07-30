##########################################################################
# Description:
    # This code applies coherent and crossing ROIs to Diffusion tensor (DT) images
    # across three samples, then outputs a histogram 
# Authors: Rhea Carlson
# Date: 12/1/2023
##########################################################################

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted
import SimpleITK as sitk

# general file path ('*' for file name that are variable across samples, F001, F002, F003)
DTimageDir='/Users/DataDir/MRI/*/DTvector.nii'
maskDir_co = '/Users/DataDir/MRI/*/mask_coherent*.nii'
maskDir_cr = '/Users/DataDir/MRI/*/mask_crossing*.nii'

resultsDir='/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/MRI/DT'

##########################################################################

# search for file paths
DT_paths=[]
for name in glob.glob(DTimageDir): 
    DT_paths.append(name)   

mask_co_paths=[]
for name in glob.glob(maskDir_co): 
    mask_co_paths.append(name)

mask_cr_paths=[]
for name in glob.glob(maskDir_cr): 
    mask_cr_paths.append(name)
    
DT_paths=natsorted(DT_paths)
mask_co_paths=natsorted(mask_co_paths)
mask_cr_paths=natsorted(mask_cr_paths)

## apply masks to DT images

def apply_DT_mask(DTpath,maskpath):  

    DTimage = sitk.ReadImage(DTpath)
    DT_arr = sitk.GetArrayFromImage(DTimage)

    mask = sitk.ReadImage(maskpath)
    mask_arr = sitk.GetArrayFromImage(mask)
    
    x = DT_arr[0] * mask_arr
    y = DT_arr[1] * mask_arr
    z = DT_arr[2] * mask_arr

    x = x[x != 0]
    y = y[y != 0]
    z = z[z != 0]

    return x,y,z

DTz_co_all=[0,0,0];DTx_co_all=[0,0,0];DTy_co_all=[0,0,0]
DTz_cr_all=[0,0,0];DTx_cr_all=[0,0,0];DTy_cr_all=[0,0,0]

for i in range(len(DT_paths)):
   [DTx_co_all[i],DTy_co_all[i],DTz_co_all[i]]=apply_DT_mask(DT_paths[i],mask_co_paths[i])
   [DTx_cr_all[i],DTy_cr_all[i],DTz_cr_all[i]]=apply_DT_mask(DT_paths[i],mask_cr_paths[i])


## Relative angle at the average angle of the ROI

## get the angle between two vectors

def unit_vector(vector):
    #Returns the unit vector of the vector
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    #Returns the angle in radians between vectors 'v1' and 'v2'::
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

# average vector (x,y,z)

avg_EV=[]
for i in range(len(DT_paths)):
    avg_EV.append([np.mean(DTx_co_all[i]),np.mean(DTy_co_all[i]),np.mean(DTz_co_all[i])])

# relative angle
ang_co=[np.zeros(len(DTx_co_all[0])),np.zeros(len(DTx_co_all[1])),np.zeros(len(DTx_co_all[2]))] # hardcoded
for i in range(len(DTx_co_all)):
    for j in range(len(DTx_co_all[i])):
        ang_co[i][j]=angle_between(avg_EV[i],(DTx_co_all[i][j],DTy_co_all[i][j], DTz_co_all[i][j]))

ang_cr=[np.zeros(len(DTx_cr_all[0])),np.zeros(len(DTx_cr_all[1])),np.zeros(len(DTx_cr_all[2]))] # hardcoded
for i in range(len(DTx_cr_all)):
    for j in range(len(DTx_cr_all[i])):
        ang_cr[i][j]=angle_between(avg_EV[i],(DTx_cr_all[i][j],DTy_cr_all[i][j], DTz_cr_all[i][j]))


#phase correction
for i in range(len(DT_paths)):
    ang_cr[i]=np.squeeze(np.where(ang_cr[i]>np.pi/2, ang_cr[i]-np.pi, ang_cr[i]))
    ang_co[i]=np.squeeze(np.where(ang_co[i]>np.pi/2, ang_co[i]-np.pi, ang_co[i]))

ang_cr=np.abs(ang_cr)
ang_co=np.abs(ang_co)

## plot histograms across 3 samples

# coherent
metrics=ang_co
metric='DT_co'
color=['c','m','y']

xtick=[0,1]
ytick=[[2,4],[1,2,3]]

fig, ax = plt.subplots(figsize=(15,20))
bin=70
min=-0
max=np.pi/2

plt.hist(x=(metrics[0]), density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[0], rwidth=0.85, label='Sample 1', color= color[0])
plt.hist(x=(metrics[1]), density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[1],rwidth=0.85, label='Sample 2', color = color[1])
plt.hist(x=(metrics[2]), density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[2],rwidth=0.85, label='Sample 3', color = color[2])

plt.title('Across Samples ' + metric, fontsize = 20, pad=50)
plt.xlim(left=min,right=max)    
plt.ylim(bottom=0,top=6)
plt.xticks(fontsize = 70)
plt.yticks(fontsize = 60)
plt.xticks(xtick) 
plt.yticks(ytick[0]) 
plt.tight_layout()

plt.tight_layout()

resultspath = os.path.join(resultsDir, 'AcrossSamples_' + metric + '.pdf')
plt.savefig(resultspath)   


#crossing
metrics=ang_cr
metric='DT_cr'
color=['c','m','y']

fig, ax = plt.subplots(figsize=(15,20))

plt.hist(x=(metrics[0]), density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[0], rwidth=0.85, label='Sample 1', color= color[0])
plt.hist(x=(metrics[1]), density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[1],rwidth=0.85, label='Sample 2', color = color[1])
plt.hist(x=(metrics[2]), density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[2],rwidth=0.85, label='Sample 3', color = color[2])

plt.title('Across Samples ' + metric, fontsize = 20, pad=50)
plt.xlim(left=min,right=max)    
plt.ylim(bottom=0,top=3.5)
plt.xticks(fontsize = 70)
plt.yticks(fontsize = 60)
plt.xticks(xtick) 
plt.yticks(ytick[1]) 

plt.tight_layout()

resultspath = os.path.join(resultsDir, 'AcrossSamples_' + metric + '.pdf')
plt.savefig(resultspath)  
