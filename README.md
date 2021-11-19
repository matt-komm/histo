# HNL analysis tools

## Installing conda (if you don't already have it)

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh
chmod +x Miniconda3-py37_4.10.3-Linux-x86_64.sh
./Miniconda3-py37_4.10.3-Linux-x86_64.sh
conda init
```

Restart the shell
Configuration files used to define samples, categories, cuts, etc. can be found in ```config```

## Installing the environment (do it only once)

NB: The ```/home/$USER``` directory storage is limited so best to clone the repository and install the conda environment on ```/vols/cms/$USER```

```
eval "$(/home/hep/$USER/miniconda3/bin/conda shell.zsh hook)"
conda create --prefix ./env python=3.9 root=6.24.6 pip -c conda-forge
conda activate ./env
pip install -r config/requirements.txt
pip install -e .
```

## Main project directories

You can set up the environment with the following. Change the first line if your conda is installed somewhere else.

```
eval "$(/home/hep/$USER/miniconda3/bin/conda shell.zsh hook)"
export HISTO_BASE_ENV=$PWD
conda activate ./env
```

Start with nanoAOD friend ntuples producing using custom [nanoAOD-tools](https://github.com/LLPDNNX/nanoAOD-tools). Additional scripts are available for single auxilliary studies with sparse documentation. Further documentation is available inside these categories:

* ```plotting```: script for signal and control region plots
* ```skim```: producing skimmed .pkl files for each signal region category to speed up subsequent studies (used for ```abcd``` and ```threshold_optimisation``` studies).
* ```abcd```: data-driven background estimation studies (work in progress)
* ```threshold_optimisation```: optimising signal region thresholds on discriminating variables
* ```limits```: Producing histograms to be used as inputs for the [CMS combine tool](https://github.com/LLPDNNX/CombineHarvester) 



