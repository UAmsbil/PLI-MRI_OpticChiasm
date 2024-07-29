############################# Description ################################

# this code finds FA and PA images for the 3 OCs, 
# applies crossing and coherent masks for the corresponding samples,
# then plots the same metric/mask for all three samples on a histogram

# inputs: FA, PA, crossing masks, coherent masks (for 3 OC samples)
# outputs: histogram image

##########################################################################

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted
#import ants <- ants was not working??
import SimpleITK as sitk

resultsDir = '/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/MRI'

#get filepaths for each FA,PA image for each sample
FAimageDir = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/blip_up_proc/*/blip_up_DMC_L0_DT_FA.nii'
PAimageDir = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/blip_up_proc/*/mapMRI/blip_up_DMC_mapmri_PA.nii'
DTimageDir='/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/blip_up_proc/*/eigen/DTvector.nii'

maskDir_co = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/blip_up_proc/*/mask_coherent*.nii'
maskDir_cr = '/xdisk/hutchinsone/rheacarlson/Data/MRI/Proc/FT_optic_chiasm_final/blip_up_proc/*/mask_crossing*.nii'

FA_paths=[]
for name in glob.glob(FAimageDir): 
    FA_paths.append(name)

PA_paths=[]
for name in glob.glob(PAimageDir): 
    PA_paths.append(name)

DT_paths=[]
for name in glob.glob(DTimageDir): 
    DT_paths.append(name)   

mask_co_paths=[]
for name in glob.glob(maskDir_co): 
    mask_co_paths.append(name)

mask_cr_paths=[]
for name in glob.glob(maskDir_cr): 
    mask_cr_paths.append(name)

FA_paths=natsorted(FA_paths)
PA_paths=natsorted(PA_paths)
DT_paths=natsorted(DT_paths)

mask_co_paths=natsorted(mask_co_paths)
mask_cr_paths=natsorted(mask_cr_paths)


# apply masks

def apply_MRI_mask(MRpath,maskpath):
    
    # # load with ants
    # MRimage = ants.image_read(MRpath)  
    # MR_arr = MRimage.numpy()
    # mask = ants.image_read(maskpath)
    # mask_arr = mask.numpy()
    # # apply mask to image
    # new_arr = MR_arr * mask_arr


    # load with simpleITK
    MRimage = sitk.ReadImage(MRpath)
    MR_arr = sitk.GetArrayFromImage(MRimage)

    mask = sitk.ReadImage(maskpath)
    mask_arr = sitk.GetArrayFromImage(mask)
    
    new_arr = MR_arr * mask_arr
    
    new_arr_flat = new_arr[new_arr != 0]

    return new_arr_flat



def apply_DT_mask(DTpath,maskpath):  

    DTimage = sitk.ReadImage(DTpath)
    DT_arr = sitk.GetArrayFromImage(DTimage)

    mask = sitk.ReadImage(maskpath)
    mask_arr = sitk.GetArrayFromImage(mask)
    
    z = DT_arr[0] * mask_arr
    x = DT_arr[1] * mask_arr
    y = DT_arr[2] * mask_arr

    z = z[z != 0]
    x = x[x != 0]
    y = y[y != 0]

    return z,x,y




FA_co_all=[0,0,0]; FA_cr_all=[0,0,0]
PA_co_all=[0,0,0]; PA_cr_all=[0,0,0]

DTz_co_all=[0,0,0];DTx_co_all=[0,0,0];DTy_co_all=[0,0,0]
DTz_cr_all=[0,0,0];DTx_cr_all=[0,0,0];DTy_cr_all=[0,0,0]

for i in range(len(FA_paths)):
   FA_co_all[i]=apply_MRI_mask(FA_paths[i], mask_co_paths[i])
   FA_cr_all[i]=apply_MRI_mask(FA_paths[i], mask_cr_paths[i])

   PA_co_all[i]=apply_MRI_mask(PA_paths[i], mask_co_paths[i])
   PA_cr_all[i]=apply_MRI_mask(PA_paths[i], mask_cr_paths[i])
   
   [DTz_co_all[i],DTx_co_all[i],DTy_co_all[i]]=apply_DT_mask(DT_paths[i],mask_co_paths[i])
   [DTz_cr_all[i],DTx_cr_all[i],DTy_cr_all[i]]=apply_DT_mask(DT_paths[i],mask_cr_paths[i])

   

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



###########################

#FA/PA AVERAGES

# FA_co_avg=(FA_co_all[0][7:1032]+FA_co_all[1]+FA_co_all[2][414:1439])/3
# FA_cr_avg=(FA_cr_all[0][61:853]+FA_cr_all[1]+FA_cr_all[2][962:1754])/3

# PA_co_avg=(PA_co_all[0][23:1016]+PA_co_all[1]+PA_co_all[2][406:1399])/3
# PA_cr_avg=(PA_cr_all[0][63:852]+PA_cr_all[1]+PA_cr_all[2][964:1753])/3


# fig, ax = plt.subplots(figsize=(15,20))
# plt.hist(x=FA_co_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=0.85, label='Sample 1', color= 'orange')
# plt.hist(x=FA_cr_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=0.85, label='Sample 2', color = 'b')
# plt.title('Avg FA', fontsize = 20, pad=30)
# plt.locator_params(axis='x', nbins=6)
# plt.yticks(fontsize = 30)
# plt.xticks(fontsize = 47)
# plt.tight_layout()
# resultspath = os.path.join(resultsDir, 'Avg_FA.pdf')
# plt.savefig(resultspath)   



# fig, ax = plt.subplots(figsize=(15,20))
# plt.hist(x=PA_co_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=0.85, label='Sample 1', color= 'orange')
# plt.hist(x=PA_cr_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=0.85, label='Sample 2', color = 'b')
# plt.title('Avg PA' , fontsize = 20, pad=30)
# plt.locator_params(axis='x', nbins=6)
# plt.yticks(fontsize = 30)
# plt.xticks(fontsize = 47)
# plt.tight_layout()
# resultspath = os.path.join(resultsDir, 'Avg_PA.pdf')
# plt.savefig(resultspath)   




###############################

########### plot DTs

#metric = 'DTx_coherent'
# metrics=DTx_co_all

# color=['c','m','y']

# fig, ax = plt.subplots(figsize=(15,20))

# bin=100
# max=1.3
# min=-1.3

# plt.hist(x=metrics[0], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor=color[0], rwidth=0.85, label='Sample 1', color= color[0])
# plt.hist(x=metrics[1], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[1],rwidth=0.85, label='Sample 2', color = color[1])
# plt.hist(x=metrics[2], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor=color[2],rwidth=0.85, label='Sample 3', color = color[2])

# plt.title('Across Samples ' + metric, fontsize = 20, pad=30)
# plt.locator_params(axis='x', nbins=4)
# plt.yticks(fontsize = 30)
# plt.xticks(fontsize = 47)

# plt.tight_layout()

# resultspath = os.path.join(resultsDir, 'AcrossSamples_' + metric + '.pdf')
# plt.savefig(resultspath)   


########## DTs AVERAGES


# FA_co_avg=(FA_co_all[0][7:1032]+FA_co_all[1]+FA_co_all[2][414:1439])/3
# FA_cr_avg=(FA_cr_all[0][61:853]+FA_cr_all[1]+FA_cr_all[2][962:1754])/3

# PA_co_avg=(PA_co_all[0][23:1016]+PA_co_all[1]+PA_co_all[2][406:1399])/3
# PA_cr_avg=(PA_cr_all[0][63:852]+PA_cr_all[1]+PA_cr_all[2][964:1753])/3


# DTx_co_avg = (DTx_co_all[0][6:1032]+DTx_co_all[1]+DTx_co_all[2][413:1439])/3
# DTy_co_avg = (DTy_co_all[0][6:1032]+DTy_co_all[1]+DTy_co_all[2][413:1439])/3
# DTz_co_avg = (DTz_co_all[0][6:1032]+DTz_co_all[1]+DTz_co_all[2][413:1439])/3


# DTx_cr_avg = (DTx_cr_all[0][61:853]+DTx_cr_all[1]+DTx_cr_all[2][962:1754])/3
# DTy_cr_avg = (DTy_cr_all[0][61:853]+DTy_cr_all[1]+DTy_cr_all[2][962:1754])/3
# DTz_cr_avg = (DTz_cr_all[0][61:853]+DTz_cr_all[1]+DTz_cr_all[2][962:1754])/3

# ########z

# fig, ax = plt.subplots(figsize=(15,20))
# bin=120
# max=1.3
# min=-1.3

# plt.hist(x=DTz_co_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=0.85, label='Sample 1', color= 'orange')
# plt.hist(x=DTz_cr_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=0.85, label='Sample 2', color = 'b')

# plt.title('Avg DTz' , fontsize = 20, pad=30)
# plt.locator_params(axis='x', nbins=6)
# plt.yticks(fontsize = 30)
# plt.xticks(fontsize = 47)

# plt.tight_layout()

# resultspath = os.path.join(resultsDir, 'Avg_DTz.pdf')
# plt.savefig(resultspath)   

# ########x


# fig, ax = plt.subplots(figsize=(15,20))

# plt.hist(x=DTx_co_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=0.85, label='Sample 1', color= 'orange')
# plt.hist(x=DTx_cr_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=0.85, label='Sample 2', color = 'b')

# plt.title('Avg DTx' , fontsize = 20, pad=30)
# plt.locator_params(axis='x', nbins=6)
# plt.yticks(fontsize = 30)
# plt.xticks(fontsize = 47)

# plt.tight_layout()

# resultspath = os.path.join(resultsDir, 'Avg_DTx.pdf')
# plt.savefig(resultspath)   

# ########y

# fig, ax = plt.subplots(figsize=(15,20))

# plt.hist(x=DTy_co_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=0.85, label='Sample 1', color= 'orange')
# plt.hist(x=DTy_cr_avg, density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=0.85, label='Sample 2', color = 'b')

# plt.title('Avg DTy' , fontsize = 20, pad=30)
# plt.locator_params(axis='x', nbins=6)
# plt.yticks(fontsize = 30)
# plt.xticks(fontsize = 47)

# plt.tight_layout()

# resultspath = os.path.join(resultsDir, 'Avg_DTy.pdf')
# plt.savefig(resultspath)  
