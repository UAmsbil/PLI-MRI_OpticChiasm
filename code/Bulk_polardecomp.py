##########################################################################
# Description: 
    # This code performs Mueller Matrix decompositiion using Lu-Chipman decomposition technique on files named for wavelengths in each image (e.g. "405", "442", etc.) 
    # for multiple images within one file
    # Outputs the "polardecomp" folder within each image file containing all Mueller matrix metrics
# Authors: Justina Bonaventura, Travis Sawyer, Rhea Carlson
# Date: 10/1/2021
##########################################################################
import sys
import os
import glob
import numpy as np
import IPython
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from pySCATMECH.mueller import *
import time

#This helps to move past poorly behaved data-
np.seterr(all='ignore')

time1 = time.time()

############################## Set Options ##############################

# This is the filepath that contains multiple PLI images
dataDir = '/Users/DataDir/PLI/20231115'

##########################################################################

# Get names of all images (subdir names)
img_name = sorted(next(os.walk(dataDir))[1])
img_name = [ x for x in img_name if 'REFERENCE' not in x and 'rep10' not in x]

print('\n\nNumber of images to analyze:', len(img_name))
print()


# Open an image
for name,num in zip(img_name,range(len(img_name))):
    folderpath = os.path.join(dataDir, name)
    print('\n(',num+1,'/',len(img_name),') Opening',folderpath)
    
    # if polardecomp folder does not already exist, perform MM decomposition
    if not os.path.exists(os.path.join(folderpath,'polardecomp')):
        
        ########################## Decomposition ########################
        
        #The code is set up to save the results to the same folder as where your mueller matrix files are.
        #This part of the code takes a long time.

        deplist=[]
        cdlist=[]
        dialist=[]
        danglist=[]
        ranglist=[]
        retlist=[]
        linretlist=[]
        optrotlist=[]
        pollist=[]

        wavelengthlist=[405, 442, 473, 543, 632]
        for i in wavelengthlist:
            wavelength=i
            wavelengthforpath="analyzed_"+str(wavelength)+"_MES"
            #This allows the code to skip over wavelengths that are not present-
            ims = glob.glob(os.path.join(folderpath, wavelengthforpath, "*.txt"))
            if not len(ims):
                print("No images are found for wavelength %s" %wavelengthforpath)
                deplist.append([0])
                cdlist.append([0])
                dialist.append([0])
                danglist.append([0])
                ranglist.append([0])
                retlist.append([0])
                linretlist.append([0])
                optrotlist.append([0])
                pollist.append([0])
                continue
            else:
                #Here the Mueller matrix text files are loaded and nan's are turned into zeros.
                print("Loading text files to image data for wavelength %s" %wavelengthforpath)
                path11=os.path.join(folderpath, wavelengthforpath, "m11.txt")
                M011=np.loadtxt(path11)
                M11 = np.where(np.isnan(M011), 0, M011)
                path12=os.path.join(folderpath, wavelengthforpath, "m12.txt")
                M012=np.loadtxt(path12)
                M12 = np.where(np.isnan(M012), 0, M012)
                path13=os.path.join(folderpath, wavelengthforpath, "m13.txt")
                M013=np.loadtxt(path13)
                M13 = np.where(np.isnan(M013), 0, M013)
                path14=os.path.join(folderpath, wavelengthforpath, "m14.txt")
                M014=np.loadtxt(path14)
                M14 = np.where(np.isnan(M014), 0, M014)
                path21=os.path.join(folderpath, wavelengthforpath, "m21.txt")
                M021=np.loadtxt(path21)
                M21 = np.where(np.isnan(M021), 0, M021)
                path22=os.path.join(folderpath, wavelengthforpath, "m22.txt")
                M022=np.loadtxt(path22)
                M22 = np.where(np.isnan(M022), 0, M022)
                path23=os.path.join(folderpath, wavelengthforpath, "m23.txt")
                M023=np.loadtxt(path23)
                M23= np.where(np.isnan(M023), 0, M023)
                path24=os.path.join(folderpath, wavelengthforpath, "m24.txt")
                M024=np.loadtxt(path24)
                M24 = np.where(np.isnan(M024), 0, M024)
                path31=os.path.join(folderpath, wavelengthforpath, "m31.txt")
                M031=np.loadtxt(path31)
                M31 = np.where(np.isnan(M031), 0, M031)
                path32=os.path.join(folderpath, wavelengthforpath, "m32.txt")
                M032=np.loadtxt(path32)
                M32 = np.where(np.isnan(M032), 0, M032)
                path33=os.path.join(folderpath, wavelengthforpath, "m33.txt")
                M033=np.loadtxt(path33)
                M33 = np.where(np.isnan(M033), 0, M033)
                path34=os.path.join(folderpath, wavelengthforpath, "m34.txt")
                M034=np.loadtxt(path34)
                M34 = np.where(np.isnan(M034), 0, M034)
                path41=os.path.join(folderpath, wavelengthforpath, "m41.txt")
                M041=np.loadtxt(path41)
                M41 = np.where(np.isnan(M041), 0, M041)
                path42=os.path.join(folderpath, wavelengthforpath, "m42.txt")
                M042=np.loadtxt(path42)
                M42 = np.where(np.isnan(M042), 0, M042)
                path43=os.path.join(folderpath, wavelengthforpath, "m43.txt")
                M043=np.loadtxt(path43)
                M43 = np.where(np.isnan(M043), 0, M043)
                path44=os.path.join(folderpath, wavelengthforpath, "m44.txt")
                M044=np.loadtxt(path44)
                M44 = np.where(np.isnan(M044), 0, M044)

                #Emply matricies are established as placeholders for data
                diatten=np.zeros(M11.shape)
                depo=np.zeros(M11.shape)
                circulardi=np.zeros(M11.shape)
                diangle=np.zeros(M11.shape)
                rang=np.zeros(M11.shape)
                ret=np.zeros(M11.shape)
                linret=np.zeros(M11.shape)
                optrot=np.zeros(M11.shape)
                pol=np.zeros(M11.shape)

                print("Performing Lu-Chipman decomposition...")
                #For each pixel in the image
                for i, row in enumerate(M11):
                    for j, col in enumerate(M11[i]):
                        #A matrix is put together from the imported text files.
                        MM=np.array([[M11[i,j], M12[i,j],M13[i,j],M14[i,j]],[M21[i,j],M22[i,j],M23[i,j],M24[i,j]],[M31[i,j],M32[i,j],M33[i,j],M34[i,j]],[M41[i,j],M42[i,j],M43[i,j],M44[i,j]]])
                        #For potential debugging-
                        #MM = np.where(MM < MM[0, 0], MM, MM[0, 0])
                        #The matrix is declared to be a Mueller matrix
                        Mule=MuellerMatrix(MM)
                        try:
                            #Parameters are able to be found using function from pySCATMech lib.
                            parametersM = CharacterizedMueller(Mule)
                            Depol=parametersM.DepolarizationCoefficient
                            Di=parametersM.Diattenuation
                            CircDi=parametersM.CircularDiattenuation
                            Diang=parametersM.DiattenuationAngle
                            ranga=parametersM.RetardanceAngle
                            reta=parametersM.Retardance
                            linreta=parametersM.LinearRetardance
                            optrota=parametersM.OpticalRotation
                            pola=parametersM.Polarizance
                        #If the matrix is not valid this will allow you to bypass and/ or debug
                        except ValueError:
                        # print(i,j)
                            Depol=0
                            Di=0
                            CircDi=0
                            Diang=0
                            ranga=0
                            reta=0
                            linreta=0
                            optrota=0
                            pola=0
                        #Calculated parameters are set to matricies.
                        depo[i,j]=Depol
                        diatten[i,j]=Di
                        circulardi[i,j]=CircDi
                        diangle[i,j]=Diang
                        rang[i,j]=ranga
                        ret[i,j]=reta
                        linret[i,j]=linreta
                        optrot[i,j]=optrota
                        pol[i,j]=pola

                #Resulting Datasets are rotated to stay consistent with other ROI code-
                deplist.append(np.rot90(np.fliplr(depo), 3))
                dialist.append(np.rot90(np.fliplr(diatten), 3))
                danglist.append(np.rot90(np.fliplr(diangle), 3))
                ranglist.append(np.rot90(np.fliplr(rang), 3))
                retlist.append(np.rot90(np.fliplr(ret), 3))
                cdlist.append(np.rot90(np.fliplr(circulardi), 3))
                linretlist.append(np.rot90(np.fliplr(linret), 3))
                optrotlist.append(np.rot90(np.fliplr(optrot), 3))
                pollist.append(np.rot90(np.fliplr(pol), 3))

        print('Saving Output')        
        # Folder called polardecomp is created within the same folder as all the mueller matrix folders-
        datapath = os.path.join(folderpath, "polardecomp") 


        if not os.path.exists(datapath):
            os.makedirs(datapath)

        #Data is saved as compressed arrays-
        for i in range(5):
            wave=wavelengthlist[i]
            np.savez_compressed(os.path.join(datapath,str(wave)), retang=ranglist[i], diatten=dialist[i], depol=deplist[i], diang=danglist[i], ret=retlist[i], cdia=cdlist[i], linret=linretlist[i], optrot=optrotlist[i], polar=pollist[i]) 

    else:
        print('polardecomp already exists!')
   

print('\nALL IMAGES FINISHED!')

# Print time to complete analysis
time2 = time.time()   
total_time = (time2 - time1)

if total_time < 60 :
    print_time = str(total_time)+' seconds'

elif total_time > 60 and total_time < 3600:
    print_time = str(total_time/60)+' minutes'

elif total_time > 3600:
    print_time = str(total_time/3600)+' hours'

total_time = (time2 - time1)/3600
print('Time to complete:', print_time)
