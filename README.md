# Contents
This repository contains all relevant code, links to the data used, masks, and expected results output from the code. PLI and dMRI data was aquired from three different OC specimens. For PLI, both thick (blockface) and thin (sectioned) OC specimens were imaged.

A representative set of MRI and PLI data used in this study is temporarily available through [Open Science Framework](https://osf.io/9pf3v/?view_only=8b66c6dc07a9496fa02caaae34919231). 

To process the PLI data, Lu-Chipman Mueller matrix decomposition was performed using the code linked [here](https://github.com/UAmsbil/PLI-MRI_OpticChiasm/blob/main/code/Bulk_polardecomp.py). PLI and MRI images were then registered via landmark registration using the code linked [here](https://github.com/UAmsbil/PLI-MRI_OpticChiasm/blob/main/code/registration_v3.mlx). Analysis was performed based on masks located [here](https://github.com/UAmsbil/PLI-MRI_OpticChiasm/tree/main/data/PLI/masks). ROI and pixelwise analysis was performed using code linked [here](https://github.com/UAmsbil/PLI-MRI_OpticChiasm/blob/main/code/histogram_across_samples.py) and [here](https://github.com/UAmsbil/PLI-MRI_OpticChiasm/blob/main/code/PLI_hist_across_samples.py).

All results from the code are [here](https://github.com/UAmsbil/PLI-MRI_OpticChiasm/tree/main/results).

# References
Lu-Chipman Mueller matrix decompsition:
1. S.-Y. Lu and R. A. Chipman, “Interpretation of Mueller matrices based on polar decomposition,” J. Opt. Soc. Am. A 13(5), 1106 (1996) [doi:10.1364/JOSAA.13.001106].
2. Python library [pySCATMECH : A Python Interface to the SCATMECH Library](https://pages.nist.gov/pySCATMECH/)

Whole brain ferret MRI:
1. E. B. Hutchinson et al., “Population based MRI and DTI templates of the adult ferret brain and tools for voxelwise analysis,” NeuroImage 152, 575–589 (2017) [doi:10.1016/j.neuroimage.2017.03.009].
