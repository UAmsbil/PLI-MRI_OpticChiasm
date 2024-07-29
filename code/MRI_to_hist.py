##########################################################################
# Description:
    # This code applies coherent and crossing ROIs to fractional anisotropy (FA) and  
    # propagator anisotropy (PA) across three OC samples, then outputs a histogram 
# Authors: Rhea Carlson
# Date: 12/1/2023
##########################################################################

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted
import SimpleITK as sitk

resultsDir = '/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/MRI'

# get filepaths for each FA,PA image for each sample
# general file path ('*' for file name that are variable across samples)
FAimageDir = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/*/blip_up_DMC_L0_DT_FA.nii'
PAimageDir = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/*/blip_up_DMC_mapmri_PA.nii'
maskDir_co = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/*/mask_coherent*.nii'
maskDir_cr = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/*/mask_crossing*.nii'

##########################################################################

FA_paths=[]
for name in glob.glob(FAimageDir): 
    FA_paths.append(name)

PA_paths=[]
for name in glob.glob(PAimageDir): 
    PA_paths.append(name)

mask_co_paths=[]
for name in glob.glob(maskDir_co): 
    mask_co_paths.append(name)

mask_cr_paths=[]
for name in glob.glob(maskDir_cr): 
    mask_cr_paths.append(name)

FA_paths=natsorted(FA_paths)
PA_paths=natsorted(PA_paths)

mask_co_paths=natsorted(mask_co_paths)
mask_cr_paths=natsorted(mask_cr_paths)


# apply masks

def apply_MRI_mask(MRpath,maskpath):
    
    # load with simpleITK
    MRimage = sitk.ReadImage(MRpath)
    MR_arr = sitk.GetArrayFromImage(MRimage)

    mask = sitk.ReadImage(maskpath)
    mask_arr = sitk.GetArrayFromImage(mask)
    
    new_arr = MR_arr * mask_arr
    
    new_arr_flat = new_arr[new_arr != 0]

    return new_arr_flat


FA_co_all=[0,0,0]; FA_cr_all=[0,0,0]
PA_co_all=[0,0,0]; PA_cr_all=[0,0,0]

for i in range(len(FA_paths)):
   FA_co_all[i]=apply_MRI_mask(FA_paths[i], mask_co_paths[i])
   FA_cr_all[i]=apply_MRI_mask(FA_paths[i], mask_cr_paths[i])

   PA_co_all[i]=apply_MRI_mask(PA_paths[i], mask_co_paths[i])
   PA_cr_all[i]=apply_MRI_mask(PA_paths[i], mask_cr_paths[i])
   

############################ plot histogram ##############################

metric = ['FA_coherent','FA_crossing']
metrics=[FA_co_all,FA_cr_all]

color=['c','m','y']

for i in range(len(metric)):

    fig, ax = plt.subplots(figsize=(15,20))

    bin=120
    max=1
    min=0

    plt.hist(x=metrics[i][0], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor=color[0], rwidth=0.85, label='Sample 1', color= color[0])
    plt.hist(x=metrics[i][1], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[1],rwidth=0.85, label='Sample 2', color = color[1])
    plt.hist(x=metrics[i][2], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[2],rwidth=0.85, label='Sample 3', color = color[2])
    plt.title('Across Samples ' + metric[i], fontsize = 20, pad=30)

    plt.yticks(fontsize = 60)
    plt.xticks(fontsize = 70)
    
    plt.xlim(left=0,right=1)    
    plt.ylim(bottom=0,top=7)
    plt.xticks([0, 0.5, 1]) 
    plt.yticks([3, 6]) 

    plt.tight_layout()

    resultspath = os.path.join(resultsDir, 'AcrossSamples_' + metric[i] + '.pdf')
    plt.savefig(resultspath)   



metrics=[PA_co_all,PA_cr_all]
metric = ['PA_coherent','PA_crossing']

for i in range(len(metric)):

    fig, ax = plt.subplots(figsize=(15,20))

    bin=120
    max=1
    min=0

    plt.hist(x=metrics[i][0], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor=color[0], rwidth=0.85, label='Sample 1', color= color[0])
    plt.hist(x=metrics[i][1], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[1],rwidth=0.85, label='Sample 2', color = color[1])
    plt.hist(x=metrics[i][2], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[2],rwidth=0.85, label='Sample 3', color = color[2])
    plt.title('Across Samples ' + metric[i], fontsize = 20, pad=30)

    plt.yticks(fontsize = 60)
    plt.xticks(fontsize = 70)
    
    plt.xlim(left=0,right=1)    
    plt.ylim(bottom=0,top=16)
    plt.xticks([0, 0.5, 1]) 
    plt.yticks([5, 10, 15]) 

    plt.tight_layout()

    resultspath = os.path.join(resultsDir, 'AcrossSamples_' + metric[i] + '.pdf')
    plt.savefig(resultspath)   

