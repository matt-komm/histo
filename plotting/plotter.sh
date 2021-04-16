#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=3:00:0
#$ -o output/$JOB_ID_$TASK_ID.out
#$ -e output/$JOB_ID_$TASK_ID.err
#$ -t 1-10


date
SGE_TASK_ID=$((SGE_TASK_ID-1))
eval "$(/vols/cms/vc1117/miniconda3/bin/conda shell.zsh hook)"
nJobs=$((10))
nVars=$((5))
yearIndex=$((SGE_TASK_ID/nJobs))
regionVarIndex=$((SGE_TASK_ID % nJobs))
years=(2016)
year=${years[yearIndex]}
regionIndex=$((regionVarIndex/nVars))
varIndex=$((regionVarIndex%nVars))
echo $regionIndex, $varIndex, $year


conda activate hnl
python -u plotter.py -r $regionIndex -v $varIndex --year $year --variable_file variables.json --region_file regions.json --output_dir plots
date
