# Histo

Classes to facilitate histogramming for the HNL analysis

* To setup the conda environment, run ```conda env create -f config/hnl.yml```
* To install the histo package, do ```pip install -e .``` (avoids problems with imports)
* An example job submission script producer is ```batch/make_hist_sub.py```
* To combine produced histograms for limits run ```./batch/merge_hists.sh```
