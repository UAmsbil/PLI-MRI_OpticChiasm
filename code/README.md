## Data Processing 
This code is run in the following order:
  1. Bulk_polardecomp.py
     * This performs Lu-chipman decoposition of the Mueller matrix
     * Outputs the polardecomp file found included with the PLI data which PLI_functions.py references in each document for accessing PLI data.
  2. registration_v3.mlx
     * This registers the middle slice from each MRI image to the corresponding PLI image for all three OC specimens for thin and thick tissue.
     * This step is necessary to perform pixelwise correlation.
     * Outputs the registered MRI images.
    
## Data Analysis 
This code is run in any order for generating histograms and pixewise correlation:
  * MRI_to_hist.py
     * Outputs histograms of MRI metrics, FA and PA, from a 3D ROI.
  * MRI_DT_to_hist.py
     * Outputs histograms of MRI metric, diffusion tensor, from a 3D ROI.
  * PLI_hist_across_samples.py
     * Outputs  histograms of PLI metrics from ROIs within the crossing and coherent regions of each OC.
     * Also outputs polarplots for angular PLI data.
  * pixelwise_corr_across_samples.py
     * This performs pixelwise Pearson correlation of dMRI and PLI data across all three specimens.
