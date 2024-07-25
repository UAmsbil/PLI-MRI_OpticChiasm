##########################################################################
# Description:
    # This code takes the histogram of ROIs within 3 coherent and 3 crossing OC PLI images
    # requires a ROI mask input, applied to all images of one type (all coherent, or all crossing)
    # It also averages values across the 3 OC specimens ROI's
# Authors: Rhea Carlson, Justina Bonaventura
# Date: 12/1/2023
##########################################################################

import numpy as np
import matplotlib.pyplot as plt
from statistics import variance
import os
import pandas as pd
from PLI_functions import *

########################### Set Filepaths #############################
MAINDIR='/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/DATA'
PLIDir = os.path.join(MAINDIR, 'data', 'PLI')
MaskDir= os.path.join(PLIDir, 'masks')

######################### Select what to run ##########################
#Run thin (0) or thick (1) slices
run_thin_or_thick=0 
#######################################################################

# get all filepaths from input
if run_thin_or_thick == 0:
    # run thin slices
    PLI_crossing = ['Crossing_Sample1_thin','Crossing_Sample2_thin', 'Crossing_Sample3_thin']
    PLI_coherent = ['Coherent_Sample1_thin', 'Coherent_Sample2_thin','Coherent_Sample3_thin']
    ROI_maskpath_crossing = os.path.join(MaskDir,'crossing_thin.tif')
    ROI_maskpath_coherent = os.path.join(MaskDir,'coherent_thin.tif')
    resultsDir = os.path.join(MAINDIR,'results','thin')
    wavelength=0 # 0 to 4 where 0=405nm and 4=632nm
    # for plotting scalar histograms
    minmax=[[0.5,1],[0,1.6],[0,0.04]] #manual
    xtick=[[0.5,0.75,1],[0, 0.8, 1.6],[0, 0.02, 0.04]]
    max_y=[40, 9, 300]
    ytick=[[18,36],[4,8],[250,125]]

elif run_thin_or_thick == 1:
    # run thick slices
    PLI_crossing = ['Crossing_Sample1_thick','Crossing_Sample2_thick', 'Crossing_Sample3_thick']
    PLI_coherent = ['Coherent_Sample1_thick', 'Coherent_Sample2_thick','Coherent_Sample3_thick']
    ROI_maskpath_crossing = os.path.join(MaskDir,'crossing_thick.tif')
    ROI_maskpath_coherent = os.path.join(MaskDir,'coherent_thick.tif')
    resultsDir = os.path.join(MAINDIR,'results','thick')
    wavelength=4 # 0 to 4 where 0=405nm and 4=632nm
    # for plotting scalar histograms
    minmax=[[0,0.5],[0,2],[0,0.06]] #manual
    xtick=[[0,0.25,0.5],[0, 1, 2],[0, 0.03, 0.06]]
    max_y=[16, 4.5, 130]
    ytick=[[7,14],[2,4],[60,120]]


############################ PLI data #################################

print('Averages across PLI metrics - ROI Mask ...')

########## Crossing ROI ##########

all_avgs_cr=[]
all_PLI_data_cr=[]
for im in PLI_crossing:
    PLIpath = os.path.join(PLIDir, im)

    # open PLI .npz folder and get data
    wvlist = getwavelist(PLIpath)
    decompDir = os.path.join(PLIpath, 'polardecomp')
    PLI_data_cr,PLI_names = listparams(decompDir, wvlist)

    ROI_mask_cr = loadmask(ROI_maskpath_crossing)
    ROI_PLI_data_cr = extract_ROI(ROI_mask_cr, PLI_data_cr) 
        
    fig, axs = plt.subplots(1, 2, sharex='col')
    plt.title(im)
    axs[0].imshow(PLI_data_cr[0][wavelength])
    axs[1].imshow(ROI_PLI_data_cr[0][wavelength])

    reshaped_cr = ROI_PLI_data_cr.reshape(9,5,1024*1024)
    all_PLI_data_cr.append(reshaped_cr)
 
metrics_list = [PLI_names.index('depolarization'),  PLI_names.index('retardance'), PLI_names.index('diattenuation'), PLI_names.index('ret_angle'), PLI_names.index('diat_angle')]

# format data we care about
all_PLI_data_cr = np.array(all_PLI_data_cr)
my_data_cr = all_PLI_data_cr[0:len(PLI_crossing),metrics_list,wavelength]

# average across images
if len(PLI_crossing)>1:
    my_avg_data_cr = np.mean(my_data_cr,0)
else:
    my_avg_data_cr = np.squeeze(my_data_cr)


########## Coherent ROI ##########

all_avgs_co=[]
all_PLI_data_co=[]
for im in PLI_coherent:
    PLIpath = os.path.join(PLIDir, im)

    # open PLI .npz folder and get data
    wvlist = getwavelist(PLIpath)
    decompDir = os.path.join(PLIpath, 'polardecomp')
    PLI_data_co,PLI_names = listparams(decompDir, wvlist)
    
    ROI_mask_co = loadmask(ROI_maskpath_coherent)
    ROI_PLI_data_co = extract_ROI(ROI_mask_co, PLI_data_co) 
    
    fig, axs = plt.subplots(1, 2, sharex='col')
    plt.title(im)
    axs[0].imshow(PLI_data_co[0][wavelength])
    axs[1].imshow(ROI_PLI_data_co[0][wavelength])
    reshaped_co = ROI_PLI_data_co.reshape(9,5,1024*1024)
    all_PLI_data_co.append(reshaped_co)
 
# format data we care about
all_PLI_data_co = np.array(all_PLI_data_co)
my_data_co = all_PLI_data_co[0:len(PLI_coherent),metrics_list,wavelength]

# average across images
if len(PLI_crossing)>1:
    my_avg_data_co = np.mean(my_data_co,0)
else:
    my_avg_data_co = np.squeeze(my_data_co)


######################### Plot Scalar Histograms ###########################

name = ['Depolarization','Retardance','Diattenuation']

for i in range(3):

    fig, ax = plt.subplots(figsize=(15,20))

    min=minmax[i][0]
    max=minmax[i][1]
    bin=200

    plt.hist(x=my_avg_data_cr[i], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=5.85, label='Crossing', color = 'b')#colors_cr[i]
    plt.hist(x=my_avg_data_co[i], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=5.85, label='Coherent', color= 'orange')#colors_co[i]

    plt.title('Histogram of Average '+ name[i]+' '+ str(wvlist[wavelength])+'nm', fontsize = 20, pad=50)
    plt.xticks(fontsize = 70)
    plt.yticks(fontsize = 60)
    plt.ylim(bottom=0,top=max_y[i])

    plt.locator_params(axis='x', nbins=3)
    plt.xticks(xtick[i]) 
    plt.yticks(ytick[i]) 
    plt.tight_layout()

    resultspath = os.path.join(resultsDir,'scalar', 'avg_' + name[i] + '.pdf')
    plt.savefig(resultspath)    


# Compare different samples across the same metric
if len(PLI_crossing)>1:

    # crossing
    for k in range(3):

        fig, ax = plt.subplots(figsize=(15,20))

        bin=200
        min=minmax[k][0]
        max=minmax[k][1]

        plt.hist(x=my_data_cr[0,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='c',rwidth=0.85, label=PLI_crossing[0], color='c')
        plt.hist(x=my_data_cr[1,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='m',rwidth=0.85, label=PLI_crossing[1], color='m')
        plt.hist(x=my_data_cr[2,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='y',rwidth=0.85, label=PLI_crossing[2], color='y')
        
        plt.title('Histogram of '+ name[k]+ ' Across Samples - Crossing', fontsize = 20,pad=30)
        plt.ylim(bottom=0,top=max_y[k])
        plt.xticks(fontsize = 70)
        plt.yticks(fontsize = 60)
        plt.xticks(xtick[k]) 
        plt.yticks(ytick[k]) 
        plt.tight_layout()

        resultspath = os.path.join(resultsDir,'scalar', 'Across_samples_crossing'+name[k]+ '.pdf')
        plt.savefig(resultspath) 

    # coherent
    for k in range(3):

        fig, ax = plt.subplots(figsize=(15,20))
        
        bin=200
        min=minmax[k][0]
        max=minmax[k][1]

        plt.hist(x=my_data_co[0,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='c', rwidth=0.85, label=PLI_coherent[0], color='c')
        plt.hist(x=my_data_co[1,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='m', rwidth=0.85, label=PLI_coherent[1], color='m')
        plt.hist(x=my_data_co[2,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='y', rwidth=0.85, label=PLI_coherent[2], color='y')

        plt.title('Histogram of '+ name[k]+ ' Across Samples - Coherent', fontsize = 20,pad=30)
        plt.ylim(bottom=0,top=max_y[k])
        plt.xticks(fontsize = 70)
        plt.yticks(fontsize = 60)
        plt.xticks(xtick[k]) 
        plt.yticks(ytick[k]) 
        plt.tight_layout()
        
        resultspath = os.path.join(resultsDir,'scalar', 'Across_samples_coherent'+name[k]+ '.pdf')
        plt.savefig(resultspath)    



############### Plot Variance of Scalar Metrics ##################

#remove nans
data_co=[[[],[],[]],[[],[],[]],[[],[],[]]]
data_cr=[[[],[],[]],[[],[],[]],[[],[],[]]]
for i in range(3):
    for j in range(3):
        data_co[i][j]= my_data_co[i][j][~np.isnan(my_data_co[i][j])]
        data_cr[i][j]= my_data_cr[i][j][~np.isnan(my_data_cr[i][j])]

#oc1: co-dep,ret,diat cr-dep,ret,diat
s1_co_var=[variance(data_co[0][0]),variance(data_co[0][1]),variance(data_co[0][2])] 
s1_cr_var=[variance(data_cr[0][0]),variance(data_cr[0][1]),variance(data_cr[0][2])] 
s2_co_var=[variance(data_co[1][0]),variance(data_co[1][1]),variance(data_co[1][2])] 
s2_cr_var=[variance(data_cr[1][0]),variance(data_cr[1][1]),variance(data_cr[1][2])] 
s3_co_var=[variance(data_co[2][0]),variance(data_co[2][1]),variance(data_co[2][2])] 
s3_cr_var=[variance(data_cr[2][0]),variance(data_cr[2][1]),variance(data_cr[2][2])] 

df=pd.DataFrame()
df['s1_co']=s1_co_var
df['s1_cr']=s1_cr_var
df['s2_co']=s2_co_var
df['s2_cr']=s2_cr_var
df['s3_co']=s3_co_var
df['s3_cr']=s3_cr_var

#df.to_csv('/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/other_PLI-OLD/thin_var.csv')

df_thin=pd.read_csv('/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/stats/thin_var.csv')
df_thick=pd.read_csv('/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/stats/thick_var.csv')

for k in range(3):

    fig = plt.subplots(figsize=(8,6))
    x = np.arange(4) 
    y1 = [df_thin['s1_co'][k],df_thick['s1_co'][k],df_thin['s1_cr'][k],df_thick['s1_cr'][k]] 
    y2 = [df_thin['s2_co'][k],df_thick['s2_co'][k],df_thin['s2_cr'][k],df_thick['s2_cr'][k]] 
    y3 = [df_thin['s3_co'][k],df_thick['s3_co'][k],df_thin['s3_cr'][k],df_thick['s3_cr'][k]]  
    width = 0.2

    plt.bar(x-0.2, y1, width,alpha=0.8,color='c',edgecolor='k') 
    plt.bar(x, y2, width,alpha=0.8,color='m',edgecolor='k') 
    plt.bar(x+0.2, y3, width,alpha=0.8,color='y',edgecolor='k') 
    plt.xticks(x, ['thin_co', 'thick_co', 'thin_cr', 'thick_cr']) 
    plt.title(name[k],pad=20)
    
    plt.locator_params(axis='y', nbins=4) 
    plt.yticks(size=35)
    ax.yaxis.offsetText.set_fontsize(24)
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    resultspath=os.path.join(resultsDir,'Variance_'+name[k]+'.pdf')
    plt.savefig(resultspath,bbox_inches="tight", dpi=300)


########################### Plot Angular Histograms  #################################

# Perform phase Correction and relative angle

name_ang = ['Retardance Angle', 'Diattenuation Angle']
metrics_list_ang=metrics_list[3:]
my_data_cr_ang = all_PLI_data_cr[0:len(PLI_crossing),metrics_list_ang,wavelength]
my_data_co_ang = all_PLI_data_co[0:len(PLI_coherent),metrics_list_ang,wavelength]

################## Crossing

#remove the nans
cr_nonans=[],[],[]
for i in range(len(PLI_crossing)): 
    for j in range(len(metrics_list_ang)):
        cr=my_data_cr_ang[i][j]
        cr_nonans[i].append(cr[np.logical_not(np.isnan(cr))])
cr_nonans=np.squeeze(cr_nonans)  

#phase correct 
for i in range(len(cr_nonans)):
    cr_nonans[i,1]=np.squeeze(np.where(cr_nonans[i,1]>0, cr_nonans[i,1]-np.pi, cr_nonans[i,1]))
    cr_nonans[i,0]=np.squeeze(np.where(cr_nonans[i,0]>np.pi/2, cr_nonans[i,0]-np.pi, cr_nonans[i,0]))

#distance from the mean    
dist_from_mean_cr=[]
for i in range(len(PLI_crossing)): 
    temp_dist_cr=[]
    for j in range(len(metrics_list_ang)):
        temp_avg_cr=np.mean(cr_nonans[i][j])
        temp_dist_cr.append(cr_nonans[i][j]-temp_avg_cr)
    dist_from_mean_cr.append(temp_dist_cr)
dist_from_mean_cr=np.squeeze(dist_from_mean_cr)
#dist_from_mean_cr=np.abs(dist_from_mean_cr)

################## Coherent

#remove the nans
co_nonans=[],[],[]
for i in range(len(PLI_coherent)): 
    for j in range(len(metrics_list_ang)):
        co=my_data_co_ang[i][j]
        co_nonans[i].append(co[np.logical_not(np.isnan(co))])
co_nonans=np.squeeze(co_nonans)  

#phase correct
for i in range(len(PLI_coherent)):
    co_nonans[i,1]=np.where(co_nonans[i,1]>0, co_nonans[i,1]-np.pi, co_nonans[i,1])
    co_nonans[i,0]=np.squeeze(np.where(co_nonans[i,0]>np.pi/2, co_nonans[i,0]-np.pi, co_nonans[i,0]))

#distancce from the mean    
dist_from_mean_co=[]
for i in range(len(PLI_coherent)): 
    temp_dist_co=[]
    for j in range(len(metrics_list_ang)):
        temp_avg_co=np.mean(co_nonans[i][j])
        temp_dist_co.append(co_nonans[i][j]-temp_avg_co)
    dist_from_mean_co.append(temp_dist_co)
dist_from_mean_co=np.squeeze(dist_from_mean_co)
#dist_from_mean_co=np.abs(dist_from_mean_co)



######################### Plot Histograms Across Samples ###########################

# Compare different samples across the same metric
if len(PLI_crossing)>1:
    
    #thin
    if run_thin_or_thick == 0:
        xtick=[-1,0,1]
        max_y=[26, 7]
        ytick=[[0,12,24],[0,3,6]]
    #thick
    if run_thin_or_thick == 1:
        xtick=[-1,0,1]
        max_y=[9, 6.2]
        ytick=[[0,4,8],[0,3,6]]

    # crossing
    for k in range(len(metrics_list_ang)):

        fig, ax = plt.subplots(figsize=(15,20))
        bin=200
        min=-np.pi/2
        max=np.pi/2

        plt.hist(x=dist_from_mean_cr[0,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='c',rwidth=0.85, label=PLI_crossing[0], color='c')
        plt.hist(x=dist_from_mean_cr[1,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='m',rwidth=0.85, label=PLI_crossing[1], color='m')
        plt.hist(x=dist_from_mean_cr[2,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='y',rwidth=0.85, label=PLI_crossing[2], color='y')
        
        plt.title('Histogram of '+ name_ang[k]+ ' Across Samples - Crossing', fontsize = 20,pad=50)
        plt.ylim(bottom=0,top=max_y[k])
        plt.xticks(fontsize = 70)
        plt.yticks(fontsize = 60)
        plt.xticks(xtick) 
        plt.yticks(ytick[k]) 
        plt.tight_layout()

        resultspath = os.path.join(resultsDir,'angular', 'Across_samples_crossing'+name_ang[k]+ '.pdf')
        plt.savefig(resultspath) 

    # coherent
    for k in range(len(metrics_list_ang)):

        fig, ax = plt.subplots(figsize=(15,20))

        plt.hist(x=dist_from_mean_co[0,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='c', rwidth=0.85, label=PLI_coherent[0], color='c')
        plt.hist(x=dist_from_mean_co[1,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='m', rwidth=0.85, label=PLI_coherent[1], color='m')
        plt.hist(x=dist_from_mean_co[2,k,:], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='y', rwidth=0.85, label=PLI_coherent[2], color='y')

        plt.title('Histogram of '+ name_ang[k]+ ' Across Samples - Coherent', fontsize = 20,pad=50)
        plt.ylim(bottom=0,top=max_y[k])
        plt.xticks(fontsize = 70)
        plt.yticks(fontsize = 60)
        plt.xticks(xtick) 
        plt.yticks(ytick[k]) 
        plt.tight_layout()
        
        resultspath = os.path.join(resultsDir, 'angular','Across_samples_coherent'+name_ang[k]+ '.pdf')
        plt.savefig(resultspath)    


# get average
my_avg_data_cr_ang=[]
my_avg_data_co_ang=[]
for j in range(len(metrics_list_ang)):
    my_avg_data_cr_ang.append((dist_from_mean_cr[0,j,:]+dist_from_mean_cr[1,j,:]+dist_from_mean_cr[2,j,:])/3)
    my_avg_data_co_ang.append((dist_from_mean_co[0,j,:]+dist_from_mean_co[1,j,:]+dist_from_mean_co[2,j,:])/3)

# plot histogram for each sample
for i in range(len(metrics_list_ang)):

    fig, ax = plt.subplots(figsize=(15,20))

    plt.hist(x=my_avg_data_cr_ang[i], density = True, bins=np.linspace(min,max,bin), alpha=0.3, edgecolor='b',rwidth=5.85, label='Crossing', color = 'b')#colors_cr[i]
    plt.hist(x=my_avg_data_co_ang[i], density = True, bins=np.linspace(min,max,bin), alpha=0.3,edgecolor='orange', rwidth=5.85, label='Coherent', color= 'orange')#colors_co[i]

    plt.title('Histogram of Average '+ name_ang[i]+' '+ str(wvlist[wavelength])+'nm', fontsize = 20, pad=50)
    plt.ylim(bottom=0,top=max_y[i])
    plt.xticks(fontsize = 70)
    plt.yticks(fontsize = 60)
    plt.xticks(xtick) 
    plt.yticks(ytick[i]) 
    plt.tight_layout()

    resultspath = os.path.join(resultsDir,'angular', 'avg_' + name_ang[i] + '.pdf')
    plt.savefig(resultspath)   


############################# Polar plots ###############################

################ Coherent:
## Retardance Angle

#Angular range
theta = np.linspace(-np.pi/2,np.pi/2, 315)

#Sort data into bins spanning the angular space-
r=np.zeros(len(theta))
r=np.histogram(my_avg_data_co_ang[0], bins=315, range=(-np.pi/2,np.pi/2))[0] #makes half of the polar plot

#Normalize
r=r/np.max(r)
thetadubs=np.linspace(-np.pi/2,3*np.pi/2, 630)
rdubs=np.concatenate((np.flip(r),np.flip(r)), axis=0) # first invert axes..

plt.figure()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(thetadubs,rdubs, color = 'b')
ax.fill(thetadubs,rdubs,color = 'b', alpha=0.25)
ax.set_rmax(1)
ax.set_rticks([0.25,0.5,0.75,1])  # Less radial ticks
plt.xticks(fontsize=12)
plt.yticks([])
ax.set_title('Retardance Angle Distribution - Coherent', va='bottom', y=1)
ax.set_theta_zero_location("N")  # theta=0 at the top - ...then rotate 90 deg
ax.set_theta_direction(-1)
plt.tight_layout()
resultspath = os.path.join(resultsDir, 'polar_plot_retang_coherent.png')
plt.savefig(resultspath,bbox_inches='tight', dpi=300) 

## Diattenuation Angle

#Angular range
theta = np.linspace(-np.pi/2,np.pi/2, 315)

#Sort data into bins spanning the angular space-
r=np.zeros(len(theta))
r=np.histogram(my_avg_data_co_ang[1], bins=315, range=(-np.pi/2,np.pi/2))[0] #makes half of the polar plot

#Normalize-
r=r/np.max(r)
thetadubs=np.linspace(-np.pi/2,3*np.pi/2, 630)
rdubs=np.concatenate((np.flip(r),np.flip(r)), axis=0) # first invert axes..

plt.figure()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(thetadubs,rdubs, color='green')
ax.fill(thetadubs,rdubs,color = 'green', alpha=0.25)
ax.set_rmax(1)
ax.set_rticks([0.25,0.5,0.75,1])  # Less radial ticks
plt.xticks(fontsize=12)
plt.yticks([])
ax.set_title('Diattenuation Angle Distribution - Coherent', va='bottom',y=1)
ax.set_theta_zero_location("N")  
ax.set_theta_direction(-1)
plt.tight_layout()
resultspath = os.path.join(resultsDir, 'polar_plot_diatang_coherent.png')
plt.savefig(resultspath,bbox_inches='tight', dpi=300) 



################ Crossing:
## Retardance Angle

#Angular range
theta = np.linspace(-np.pi/2,np.pi/2, 315)

#Sort data into bins spanning the angular space-
r=np.zeros(len(theta))
r=np.histogram(my_avg_data_cr_ang[0], bins=315, range=(-np.pi/2,np.pi/2))[0] #makes half of the polar plot

#Normalize
r=r/np.max(r)
thetadubs=np.linspace(-np.pi/2,3*np.pi/2, 630)
rdubs=np.concatenate((np.flip(r),np.flip(r)), axis=0) # first invert axes..

plt.figure()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(thetadubs,rdubs, color = 'b')
ax.fill(thetadubs,rdubs,color = 'b', alpha=0.25)
ax.set_rmax(1)
ax.set_rticks([0.25,0.5,0.75,1])  # Less radial ticks
plt.xticks(fontsize=12)
plt.yticks([])
ax.set_title('Retardance Angle Distribution - Crossing', va='bottom', y=1)
ax.set_theta_zero_location("N")  # theta=0 at the top - ...then rotate 90 deg
ax.set_theta_direction(-1)
plt.tight_layout()
resultspath = os.path.join(resultsDir, 'polar_plot_retang_crossing.png')
plt.savefig(resultspath,bbox_inches='tight', dpi=300) 


## Diattenuation Angle

#Angular range
theta = np.linspace(-np.pi/2,np.pi/2, 315)

#Sort data into bins spanning the angular space-
r=np.zeros(len(theta))
r=np.histogram(my_avg_data_cr_ang[1], bins=315, range=(-np.pi/2,np.pi/2))[0] #makes half of the polar plot

#Normalize-
r=r/np.max(r)
thetadubs=np.linspace(-np.pi/2,3*np.pi/2, 630)
rdubs=np.concatenate((np.flip(r),np.flip(r)), axis=0) # first invert axes..

plt.figure()

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(thetadubs,rdubs, color='green')
ax.fill(thetadubs,rdubs,color = 'green', alpha=0.25)
ax.set_rmax(1)
ax.set_rticks([0.25,0.5,0.75,1])  # Less radial ticks
plt.xticks(fontsize=12)
plt.yticks([])
ax.set_title('Diattenuation Angle Distribution - Crossing', va='bottom',y=1)
ax.set_theta_zero_location("N")  # theta=0 at the top - ...then rotate 90 deg
ax.set_theta_direction(-1)
plt.tight_layout()
resultspath = os.path.join(resultsDir, 'polar_plot_diatang_crossing.png')
plt.savefig(resultspath,bbox_inches='tight', dpi=300) 