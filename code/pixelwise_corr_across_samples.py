##########################################################################
# Description:
    # This code is run after Matlab script that registers MRI to PLI
    # It performs pixelwise correlation between select pre-registered MRI .png images and all named PLI data for one image
    # inputs: pre-registered MRI images, PLI polardecomp file (only crossing), PLI area mask
# Authors: Rhea Carlson, Justina Bonaventura
# Date: 12/1/2023
##########################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image
import os
import pandas as pd
import seaborn as sns
from PLI_functions import *

########################### Set Filepaths ###############################

MAINDIR='/xdisk/hutchinsone/rheacarlson/Workspace/PLI_hist/DATA'
PLIDir = os.path.join(MAINDIR, 'data', 'PLI')
MaskDir= os.path.join(PLIDir, 'masks')

######################### Select what to run ##########################
#Run thin (0) or thick (1) slices
run_thin_or_thick=1
#######################################################################

# get all filepaths from input
if run_thin_or_thick == 0:
    # run thin slices
    state="thin"
    MRI_slice_Dir=os.path.join(MAINDIR,"results","registration","thin")
    PLI_images=['Crossing_Sample1_thin','Crossing_Sample2_thin','Crossing_Sample3_thin']
    ROI_maskpath=os.path.join(PLIDir,"masks","whole_thin.tif")
    subset1=['FA','PA','DTvector0','DTvector1','DTvector2','rel_DTvector','depolarization_405', 'retardance_405','diattenuation_405', 'ret_angle_405','diat_angle_405','rel_ret_angle_405','rel_diat_angle_405']
    subset2=['FA','PA','depolarization_405', 'retardance_405','diattenuation_405']
    subset3=['rel_DTvector', 'rel_ret_angle_405','rel_diat_angle_405']
    resultspath=os.path.join(MAINDIR,"results","correlation","thin")
elif run_thin_or_thick == 1:
    # run thick slices
    state="thick"
    MRI_slice_Dir=os.path.join(MAINDIR,"results","registration","thick")
    PLI_images=['Crossing_Sample1_thick','Crossing_Sample2_thick','Crossing_Sample3_thick']
    ROI_maskpath=os.path.join(PLIDir,"masks","whole_thick.tif")
    subset1=['FA','PA','DTvector0','DTvector1','DTvector2','rel_DTvector','depolarization_632', 'retardance_632','diattenuation_632', 'ret_angle_632','diat_angle_632','rel_ret_angle_632','rel_diat_angle_632']
    subset2=['FA','PA','depolarization_632', 'retardance_632','diattenuation_632']
    subset3=['rel_DTvector', 'rel_ret_angle_632','rel_diat_angle_632']
    resultspath=os.path.join(MAINDIR,"results","correlation","thick")

map = 'RdBu_r' # colormap used for corr plots

############################ PLI data #################################

a=0
for im in PLI_images:
    sample_PLI_data=[]
    PLIpath = os.path.join(PLIDir, im)
    #print(PLIpath)

    # open PLI .npz folder and get data
    wvlist = getwavelist(PLIpath)
    decompDir = os.path.join(PLIpath, 'polardecomp')
    PLI_data,PLI_types = listparams(decompDir, wvlist)

    ROI_mask = loadmask(ROI_maskpath)
    ROI_PLI_data = extract_ROI(ROI_mask, PLI_data) 

    sample_PLI_data.append(ROI_PLI_data)
    sample_PLI_data=np.squeeze(sample_PLI_data)

    ############################ MRI data ##################################

    # filter/order relevant MRI images - could just use glob...
    MRI_slice_paths = os.listdir(MRI_slice_Dir)
    MRI_slice_paths = [x for x in MRI_slice_paths if '.png' in x and not 'error' in x and not 'avg-z' in x and im in x]
    MRI_slice_paths=sorted(MRI_slice_paths)

    MRI_data=[]
    MRI_type=[]
    for path in MRI_slice_paths:
        # load MR image from .png
        fullpath = os.path.join(MRI_slice_Dir,path)
        print(fullpath)
        temp_image = image.imread(fullpath)

        if path == MRI_slice_paths[3]:
            my_image=temp_image
        
        # apply mask to MR image 
        temp_image = temp_image * ROI_mask

        MRI_data.append(temp_image)
        MRI_type.append(os.path.splitext(path)[0])
    
    plt.subplot(2,2,1)
    plt.imshow(PLI_data[0][0])
    plt.title(im, pad=10)  
    plt.subplot(2,2,2)
    plt.imshow(sample_PLI_data[0][0])
    plt.subplot(2,2,3)
    plt.imshow(my_image, cmap='gray')
    plt.subplot(2,2,4)
    plt.imshow(MRI_data[3], cmap='gray')
    plt.show()

    ######################### Create df ################################

    # flatten arrays to prep for df creation
    sample_PLI_data = sample_PLI_data.reshape((9,5,1024*1024))
    
    MRI_data_flat=MRI_data.copy()
    for i in range(len(MRI_data)):
        MRI_data_flat[i] = MRI_data[i].ravel()

    # create df
    df = pd.DataFrame()
    for i in range(len(MRI_type)):
        df[MRI_type[i].split('_')[0]] = MRI_data_flat[i]

    for i in range(len(PLI_types)):
        for j in range(len(wvlist)):
            df[PLI_types[i]+'_'+str(wvlist[j])] = sample_PLI_data[i][j]

    #remove rows containing NaNs
    df_noZeros = df.copy().dropna()

    #remove rows where PA or FA = 0
    MRI_zeros = df[(df['FA'] == 0) | (df['PA'] == 0)].index
    df_noZeros.drop(MRI_zeros, inplace=True)
    df_noZeros.reset_index()

    ######### get relative eigenvector from DT x,y,z 

    def unit_vector(vector):
        #Returns the unit vector of the vector
        return vector / np.linalg.norm(vector)

    def angle_between(v1, v2):
        #Returns the angle in radians between vectors 'v1' and 'v2'::
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    # get average eigenvector for each DT x,y,z
    avg_EV=[]

    avg_EV.append([np.mean(df_noZeros.loc[:,'DTvector0']),np.mean(df_noZeros.loc[:,'DTvector1']),np.mean(df_noZeros.loc[:,'DTvector2'])])
    avg_EV=np.squeeze(avg_EV)

    # get the relative angle between average DT and DT xyz
    relative_ang=[]    
    DTvector0=df_noZeros.loc[:,'DTvector0'].to_numpy()
    DTvector1=df_noZeros.loc[:,'DTvector1'].to_numpy()
    DTvector2=df_noZeros.loc[:,'DTvector2'].to_numpy()

    for i in range(len(DTvector0)):
        relative_ang.append(angle_between(avg_EV,[DTvector0[i],DTvector1[i],DTvector2[i]]))

    df_noZeros['rel_DTvector']=relative_ang

    ########## get relative diattenuation and retardance!!!

    retang=df_noZeros.loc[:,subset1[-4]] # positions are hardcoded
    diatang=df_noZeros.loc[:,subset1[-3]]

    #phase correct
    retang=np.squeeze(np.where(retang>np.pi/2, retang-np.pi, retang))
    diatang=np.squeeze(np.where(diatang>0, diatang-np.pi, diatang))

    #distancce from the mean    
    avg_retang=np.mean(retang)
    avg_diatang=np.mean(diatang)
    
    dist_retang=np.squeeze(retang-avg_retang)
    dist_diatang=np.squeeze(diatang-avg_diatang)

    df_noZeros['rel_'+subset1[-4]]=dist_retang 
    df_noZeros['rel_'+subset1[-3]]=dist_diatang

    # save dfs
    if a==0:
        df1=df_noZeros
    if a==1:
        df2=df_noZeros
    if a==2:
        df3=df_noZeros

    a=a+1


###################### get corrs for each sample & average ######################

corr1 = df1.corr(method='pearson')
corr1['DTvector2'].fillna(0, inplace=True)
corr1.loc['DTvector2'].fillna(0, inplace=True)
corr2 = df2.corr(method='pearson')
corr2['DTvector2'].fillna(0, inplace=True)
corr2.loc['DTvector2'].fillna(0, inplace=True)
corr3 = df3.corr(method='pearson')
corrM = pd.concat([corr1,corr2,corr3]).groupby(level=0).mean() # average df values

####################### Plot Single Sample Corrs #############################

all_subsets=[subset1, subset2, subset3]

def single_image_corr(corr,subset, state, name, resultspath, all_subsets):

    if subset==all_subsets[1] or subset==all_subsets[2]: 
        if 'FA' in subset[0]:
            sz=34
            resultspath=os.path.join(resultspath,'scalar')
            char='scalar'
        elif 'DT' in subset[0]:
            sz=44
            resultspath=os.path.join(resultspath,'angular')
            char='angular'
    elif subset==all_subsets[0]: 
        sz=20
        resultspath=os.path.join(resultspath,'both')
        char='both'

    else:
        sz=20
        resultspath=resultspath
        char='uncategorized'
        
    new_corr = corr[subset].copy()
    new_corr=new_corr.loc[subset]
    plt.figure(figsize=[15,12])
    ax = sns.heatmap(new_corr, annot = True, fmt='.2f',linewidths = .1, cmap=map, 
                    cbar_kws={'label': 'Pearson Correlation Coeff.', 'shrink': 0.9}, annot_kws={"size": sz},  vmin = -1, vmax = 1)
    ax.figure.axes[-1].yaxis.label.set_size(15) #changes the size of the colorbar label
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.tick_params(axis='x', rotation=90)
    plt.title(name+'_'+state+'_'+char, fontsize = 20, pad=20)
    plt.tight_layout()


    plt.savefig(os.path.join(resultspath, name+'_'+state+'_selection_pixelwise.pdf'))

    return


# chose what subset to plot

# all 
single_image_corr(corr1,subset1, state,'F001', resultspath,all_subsets) 
single_image_corr(corr2,subset1, state,'F002', resultspath,all_subsets) 
single_image_corr(corr3,subset1, state,'F003', resultspath,all_subsets) 
single_image_corr(corrM,subset1, state,'AVG', resultspath,all_subsets) 

# scalar
single_image_corr(corr1,subset2, state,'F001', resultspath,all_subsets) 
single_image_corr(corr2,subset2, state,'F002', resultspath,all_subsets) 
single_image_corr(corr3,subset2, state,'F003', resultspath,all_subsets) 
single_image_corr(corrM,subset2, state,'AVG', resultspath,all_subsets) 

# angular
single_image_corr(corr1,subset3, state,'F001', resultspath,all_subsets) 
single_image_corr(corr2,subset3, state,'F002', resultspath,all_subsets) 
single_image_corr(corr3,subset3, state,'F003', resultspath,all_subsets) 
single_image_corr(corrM,subset3, state,'AVG', resultspath,all_subsets) 



########################## Plot Averages - ALL ################################

# Plot correlation matrix
plt.figure(figsize=[20,10])
ax = sns.heatmap(corrM, annot = True, fmt='.2f',linewidths = .1, cmap=map, 
                cbar_kws={'label': 'Pearson Correlation Coeff.', 'shrink': 0.8}, annot_kws={"size": 6},  vmin = -1, vmax = 1)
ax.figure.axes[-1].yaxis.label.set_size(15) #changes the size of the colorbar label
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax.tick_params(axis='x', rotation=90)
plt.title('Average_'+state, fontsize = 20, pad=20)
plt.tight_layout()

# check if resultspath exists
if not os.path.isdir(resultspath):
    print('Creating results path...')
    os.makedirs(resultspath)

# save out plot
plt.savefig(os.path.join(resultspath, 'Avg_'+state+'_ALLmetrics_pixelwise.pdf'))



########################## Scatter plots ################################

# plot subset2 and subset3 correlations on scatterplots

for a in [0,1]:
    for b in [3,4]:

        plt.figure(figsize=[5,5])
        plt.scatter(df1[subset2[a]], df1[subset2[b]],s=0.01,c='c',alpha=0.3)
        plt.scatter(df2[subset2[a]], df2[subset2[b]],s=0.01,c='m',alpha=0.3)
        plt.scatter(df3[subset2[a]], df3[subset2[b]],s=0.01,c='y',alpha=0.3)

        plt.xlabel(subset2[a])
        plt.ylabel(subset2[b])

        pearson_corr1 = df1[subset2[a]].corr(df1[subset2[b]])
        pearson_corr2 = df2[subset2[a]].corr(df2[subset2[b]])
        pearson_corr3 = df3[subset2[a]].corr(df3[subset2[b]])

        l=0.4
        corr_text1 = "R={:.2f}".format(pearson_corr1)
        plt.annotate(corr_text1, xy=(1.05, 0.15+l), xycoords='axes fraction',color='c')
        corr_text2 = "R={:.2f}".format(pearson_corr2)
        plt.annotate(corr_text2, xy=(1.05, 0.1+l), xycoords='axes fraction',color='m')
        corr_text3 = "R={:.2f}".format(pearson_corr3)
        plt.annotate(corr_text3, xy=(1.05, 0.05+l), xycoords='axes fraction',color='y')

        plt.savefig(os.path.join(resultspath, 'scatter', subset2[a]+'-'+subset2[b]+'.png'))

for b in [1,2]:

    plt.figure(figsize=[5,5])
    plt.scatter(df1[subset3[0]], df1[subset3[b]],s=0.01,c='c',alpha=0.3)
    plt.scatter(df2[subset3[0]], df2[subset3[b]],s=0.01,c='m',alpha=0.3)
    plt.scatter(df3[subset3[0]], df3[subset3[b]],s=0.01,c='y',alpha=0.3)

    plt.xlabel(subset3[0])
    plt.ylabel(subset3[b])

    pearson_corr1 = df1[subset3[0]].corr(df1[subset3[b]])
    pearson_corr2 = df2[subset3[0]].corr(df2[subset3[b]])
    pearson_corr3 = df3[subset3[0]].corr(df3[subset3[b]])

    l=0.4
    corr_text1 = "R={:.2f}".format(pearson_corr1)
    plt.annotate(corr_text1, xy=(1.05, 0.15+l), xycoords='axes fraction',color='c')
    corr_text2 = "R={:.2f}".format(pearson_corr2)
    plt.annotate(corr_text2, xy=(1.05, 0.1+l), xycoords='axes fraction',color='m')
    corr_text3 = "R={:.2f}".format(pearson_corr3)
    plt.annotate(corr_text3, xy=(1.05, 0.05+l), xycoords='axes fraction',color='y')

    plt.savefig(os.path.join(resultspath, 'scatter', subset3[0]+'-'+subset3[b]+'.png'))
