# Contents
This repository contains all relevant code, links to the data used, masks, and expected results output from the code. PLI and dMRI data was aquired from three different OC specimens. For PLI, both thick (blockface) and thin (sectioned) OC specimens were imaged.

To process the PLI data, Lu-Chipman Mueller matrix decomposition was performed using the code linked [here](https://github.com/rcarl0/OC_PLI-dMRI/blob/main/code/Bulk_polardecomp.py). PLI and MRI images were then registered via landmark registration using the code linked [here](https://github.com/rcarl0/OC_PLI-dMRI/blob/main/code/registration_v3.mlx). Region of interest (ROI) analysis was performed based on masks located [here](https://github.com/rcarl0/OC_PLI-dMRI/tree/main/data/PLI/masks). ROI and pixelwise analysis was performed using code linked [here](https://github.com/rcarl0/OC_PLI-dMRI/blob/main/code/histogram_across_samples.py) and [here](https://github.com/rcarl0/OC_PLI-dMRI/blob/main/code/pixelwise_corr_across_samples.py).

All PLI and MRI data can be accessed through the following DOI:

# References
Lu-Chipman Mueller matrix decompsition:
1. S.-Y. Lu and R. A. Chipman, “Interpretation of Mueller matrices based on polar decomposition,” J. Opt. Soc. Am. A 13(5), 1106 (1996) [doi:10.1364/JOSAA.13.001106].
2. Python library [pySCATMECH : A Python Interface to the SCATMECH Library](https://pages.nist.gov/pySCATMECH/)
