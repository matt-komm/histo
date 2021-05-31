# Tools to make histograms for HNL analysis and more:

## Classes to facilitate histogramming for the HNL analysis: 

* To setup the conda environment, run ```conda env create -f config/hnl.yml```
* To install the histo package, do ```pip install -e .``` (avoids problems with imports)
* Configuration files used to define samples, categories, cuts, etc. are stored in ```config```

## Studies this is used for:
* Producing histograms to be used as inputs for the CMS combine tool ```limits```
* General plotting script (SR/CR) ```plotting```
* Scripts for ABCD studies (background estimation) ```abcd```
* Making HNL kinematic plots ```hnl_kinematics```
* Performing lepton efficiency studies ```lepton_efficiency```
* Plotting systematic uncertainty profiles ```uncertainty_profiles```