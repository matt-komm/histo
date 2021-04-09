
categories_2l=[
    "mumu_OS_displaced", 
    "mumu_SS_displaced",
    "ee_OS_displaced",
    "ee_SS_displaced",
    "mue_OS_displaced",
    "mue_SS_displaced",
    "emu_OS_displaced",
    "emu_SS_displaced",
    "mumu_OS_prompt",
    "mumu_SS_prompt",
    "ee_OS_prompt",
    "ee_SS_prompt",
    "mue_OS_prompt",
    "mue_SS_prompt",
    "emu_OS_prompt",
    "emu_SS_prompt",
    ]

categories_1l = ["mu_single", "e_single"]

ncats_2l = len(categories_2l)
ncats_1l = len(categories_1l)


regions = ["A", "B", "C", "D"]
nregions = len(regions)

procsbkg = []
procsHNL = []
procsData = []
years = ["2016", "2017", "2018"]
leptons = ["1", "2"]

import yaml
with open("samples.yml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)

for l in samples_dict:
    if "HNL" in l:
        procsHNL.append(l)
    elif l == "muon" or l == "electron": 
        continue
        procsData.append(l)
    else:
        continue
        procsbkg.append(l)

nprocsbkg = len(procsbkg)
nprocsHNL = len(procsHNL)
nprocsData = len(procsData)

njobs_2l = nregions*nprocsbkg*ncats_2l + nprocsHNL*ncats_2l + int(nregions*ncats_2l*nprocsData/2)
print(f"The number of 2l jobs is {njobs_2l}.")

njobs_1l = nregions*nprocsbkg*ncats_1l + nprocsHNL*ncats_1l + int(nregions*ncats_1l*nprocsData/2)
print(f"The number of 1l jobs is {njobs_1l}.")

for lepton in leptons:
    if lepton == "1":
        categories = categories_1l
        njobs = njobs_1l
    else:
        categories = categories_2l
        njobs = njobs_2l
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
        with open(f"subHists{lepton}_{year}.sh", "w") as f:
            f.write(sub_string)
            f.write("date\n")
            f.write('eval "$(/vols/cms/$USER/miniconda3/bin/conda shell.bash hook)"; conda activate hnl\n')
            f.write("JOBS=(")
            for proc in procsbkg:
                for region in regions:
                    for category in categories:
                        f.write(f'"python -u make_hists.py --proc {proc} --category {category} --region {region} --year {year} --leptons {lepton}"\n')
            
            for proc in procsHNL:
                for category in categories:
                    f.write(f'"python -u make_hists.py --proc {proc} --category {category} --region D --year {year} --leptons {lepton}"\n')

            for proc in procsData:
                for region in regions:
                    for category in categories:
                        if proc == "muon":
                            if category.startswith("mu"):
                                f.write(f'"python -u make_hists.py --data --proc {proc} --category {category} --region {region} --year {year} --leptons {lepton}"\n')
                        elif proc == "electron":
                            if category.startswith("e"):
                                f.write(f'"python -u make_hists.py --data --proc {proc} --category {category} --region {region} --year {year} --leptons {lepton}"\n')

            f.write(")\n")
            f.write("echo ${JOBS[$SGE_TASK_ID-1]}\n")

            f.write("${JOBS[$SGE_TASK_ID-1]}\n")
            f.write("date\n")
