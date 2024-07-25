##########################################################################
# Description: 
    # The purpose of this code is to define functions for: 
    # Loading PLI data from polardecomp file
    # Loading mask for PLI data
    # Applying mask to PLI data
# Authors: Justina Bonaventura, Rhea Carlson
# Date: 10/1/2021   
##########################################################################

import numpy as np
import os
from PIL import Image
import cv2
import glob
import matplotlib.pyplot as plt

##########################################################################

#This function checks the folder different wavelength data and returns a list of what is present- 
def getwavelist(folderpath):
    wavelist=[]
    waveopts=[405, 442, 473, 543, 632]
    for wv in waveopts:
        wavelengthforpath="analyzed_"+str(wv)+"_MES"
        ims = glob.glob(os.path.join(folderpath, wavelengthforpath, "*.txt"))
        if len(ims):
            wavelist.append(wv)
    return(wavelist)


# This function loads data from compressed files and organizes it all into a list of lists-
def listparams(folderpath, wavelist):
    waveopts=[405, 442, 473, 543, 632]
    
    filler= np.array([0, 0, 0, 0])
    depolarization=[]
    circ_diatten=[]
    diattenuation=[]
    diat_angle=[]
    ret_angle=[]
    retardance=[]
    linear_ret=[]
    optical_rotation=[]
    polarizance=[]

    for i in range(len(waveopts)):
        if waveopts[i] in wavelist:  
            # images are flipped and rotated upon loading
            wave=waveopts[i]
            wavefile=np.load(os.path.join(folderpath, str(wave)+".npz"))
            ret_angle.append(np.fliplr(np.rot90(wavefile['retang'],1)))
            diattenuation.append(np.fliplr(np.rot90(wavefile['diatten'],1)))
            depolarization.append(np.fliplr(np.rot90(wavefile['depol'],1)))
            diat_angle.append(np.fliplr(np.rot90(wavefile['diang'],1)))
            retardance.append(np.fliplr(np.rot90(wavefile['ret'],1)))           
            circ_diatten.append(np.fliplr(np.rot90(wavefile['cdia'],1)))
            linear_ret.append(np.fliplr(np.rot90(wavefile['linret'],1)))
            optical_rotation.append(np.fliplr(np.rot90(wavefile['optrot'],1)))
            polarizance.append(np.fliplr(np.rot90(wavefile['polar'],1)))   

        else: 
            diattenuation.append(filler)
            depolarization.append(filler)
            circ_diatten.append(filler)
            diat_angle.append(filler)
            ret_angle.append(filler) 
            retardance.append(filler) 
            linear_ret.append(filler) 
            optical_rotation.append(filler) 
            polarizance.append(filler)
    
    # organized into the nine different PLI metrics, subdivided by 5 wavelengths for each metrics 
    listlist=[depolarization, circ_diatten, diattenuation, diat_angle, ret_angle, retardance, linear_ret, optical_rotation, polarizance]
    feature_names =['depolarization', 'circ_diatten', 'diattenuation', 'diat_angle', 'ret_angle', 'retardance', 'linear_ret', 'optical_rotation', 'polarizance']
    
    return(listlist, feature_names)


#This function loads a mask and makes sure it is properly scaled 
def loadmask(maskfolderpath): 
    mask = np.asarray(Image.open(maskfolderpath))
    #Standardize Mask
    mask=mask/np.max(mask)-np.min(mask)
    if mask.shape[0]>mask.shape[1]:
        sqmask=mask[:mask.shape[1], :]
    else:
        sqmask=mask[:,:mask.shape[0]]
    mask=np.abs(cv2.resize(sqmask,(1024,1024))-1)
    mask=np.where(mask==1,mask,np.nan)

    return(mask)


# This function applies a mask to all data
def extract_ROI(mask, listlist):
   
    waveopts=[405, 442, 473, 543, 632]
    ROI_data=np.full_like(listlist,0)

    for i in range(len(listlist)):
        for j in range(len(listlist[i])):
            temp_data=[]

            if len(listlist[i][j])>5:
                temp_data = listlist[i][j]*mask
            else:
                temp_data = np.zeros((1024,1024))

            ROI_data[i][j] = temp_data
                
    return(ROI_data)
