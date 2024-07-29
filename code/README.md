## Data Processing
This code is run in the following order:
  1. Bulk_polardecomp.py
     * This outputs the polardecomp file found included with the PLI data.
     * This file is referenced by PLI_functions.py in each document
  2. registration_v3.mlx
     * This registers the middle slice from the MRI image to the corresponding PLI image for all three OC specimens for thin and thick tissue.
     * Outputs the registered MRI images.
    
## Data Analysis
This code is run in any order:
  * histogram_across_samples.py
     * Outputs  histograms of PLI metrics from ROIs within the crossing and coherent regions of each OC
     * Outputs polarplots for angular PLI data 
  * pixelwise_corr_across_samples.py
     * Pearson correlation of dMRI and PLI data across all three specimens
     * Also 


  
