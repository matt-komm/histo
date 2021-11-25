
categories=[
    "mumu_OS", 
    "mumu_SS",
    "ee_OS",
    "ee_SS",
    "mue_OS",
    "mue_SS",
    "emu_OS",
    "emu_SS",
    ]

ncats = len(categories)

regions = ["A", "B", "C", "D"]
nregions = len(regions)

procsbkg = []
procsHNL = []
procsData = []
years = ["2016", "2017", "2018"]

import yaml
with open("../config/samples.yml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)

for l in samples_dict:
    if "HNL" in l:
        procsHNL.append(l)
    elif l == "muon" or l == "electron": 
        procsData.append(l)
    else:
        continue
        procsbkg.append(l)

nprocsbkg = len(procsbkg)
nprocsHNL = len(procsHNL)
nprocsData = len(procsData)

njobs = nregions*(nprocsbkg+nprocsHNL)*ncats + int(nregions*ncats*nprocsData/2)
print(f"The number of jobs is {njobs}.")

for year in years:
    sub_string = """#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=3:00:00
#$ -l h_vmem=24G
#$ -o output/$JOB_ID_$TASK_ID.out
#$ -e output/$JOB_ID_$TASK_ID.err
#$ -t 1-{}
""".format(njobs)
    with open(f"subHists_{year}.sh", "w") as f:
        f.write(sub_string)
        f.write("date\n")
        f.write('eval "$(/vols/cms/$USER/miniconda3/bin/conda shell.bash hook)"; conda activate hnl\n')
        f.write("JOBS=(")
        for proc in procsbkg+procsHNL:
            for region in regions:
                for category in categories:
                    f.write(f'"python -u make_hists.py --proc {proc} --category {category} --region {region} --year {year} "\n')

        for proc in procsData:
            for region in regions:
                for category in categories:
                    if proc == "muon" and category.startswith("mu"):
                        f.write(f'"python -u make_hists.py --data --proc {proc} --category {category} --region {region} --year {year} "\n')
                    elif proc == "electron" and category.startswith("e"):
                            f.write(f'"python -u make_hists.py --data --proc {proc} --category {category} --region {region} --year {year} "\n')

        f.write(")\n")
        f.write("echo ${JOBS[$SGE_TASK_ID-1]}\n")

        f.write("${JOBS[$SGE_TASK_ID-1]}\n")
        f.write("date\n")
