# Make limit histograms 

* The main drive is ```make_hists.py```. Start with processed [nanoAOD-tools](https://github.com/LLPDNNX/nanoAOD-tools) files.
* An SGE batch job submission script producer is ```make_hist_sub.py```
* To combine produced histograms for limits run ```merge_hists.sh```
* To run limits consult [CombineHarvester](https://github.com/LLPDNNX/CombineHarvester) which creates cards to run with combine
