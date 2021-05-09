# Tools to make histograms for HNL analysis and more

## Classes to facilitate histogramming for the HNL analysis

* To setup the conda environment, run ```conda env create -f config/hnl.yml```
* To install the histo package, do ```pip install -e .``` (avoids problems with imports)

## Producing histograms to be used as inputs for the CMS combine tool
* The main drive is ```limits/make_hists.py```. Start with processed nanoAOD-tools (Friend) files.
* An example job submission script producer is ```limits/make_hist_sub.py```
* To combine produced histograms for limits run ```./limits/merge_hists.sh```
* To make yields table and plot look at ```limits/yields.py```.

## Some basic HNL kinematic plots
* Have a look at hnl_kinematics directory

## Lepton efficiency studies
* Have a look at lepton_efficiency category