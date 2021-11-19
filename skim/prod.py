import os
import json
import yaml
import argparse
import numpy as np
import pandas as pd
from histo import Sample, Process

parser = argparse.ArgumentParser()

parser.add_argument("--category", action="store", default="mumu_OS_displaced")
parser.add_argument("--output_dir", action="store", default="/vols/cms/$USER/AN-19-207/friend_skim/")
parser.add_argument("--year", action="store", default="2016")
parser.add_argument("--test", action="store_true", dest="one_file", default=False)

args = parser.parse_args()
category_name = args.category

categories_file = os.path.expandvars("$HISTO_BASE_PATH/config/categories_2l_inclusive.json")
print(categories_file)
ntuple_path = f"/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21"
output_dir = os.path.expandvars(f"{args.output_dir}")
one_file = args.one_file
year = args.year

lumi = {"2016": 35.92, "2017": 41.53, "2018": 59.74}
with open("../config/samples.yml", 'r') as samples_file:
    samples_dict = yaml.safe_load(samples_file)

with open(categories_file, 'r') as fp:
    weight_dict = json.load(fp)

processes = []

text_dict = {
    "HNL_majorana_pt20_ctau1p0e-03_massHNL16p0_Vall1p097e-02": "HNL",
    "HNL_majorana_pt20_ctau1p0e-01_massHNL12p0_Vall2p314e-03": "HNL",
    "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": "HNL",
    "HNL_majorana_pt20_ctau1p0e01_massHNL8p0_Vall6p702e-04": "HNL",
    "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": "HNL",
    "HNL_majorana_pt20_ctau1p0e02_massHNL6p0_Vall4p597e-04": "HNL",
    "HNL_majorana_pt20_ctau1p0e03_massHNL2p0_Vall2p871e-03": "HNL",
    "HNL_majorana_pt20_ctau1p0e03_massHNL3p0_Vall9p825e-04": "HNL",
    "wjets": "W+Jets",
    "dyjets": "DY+Jets",
    "topbkg": "t#bar{t}/st",
    "vgamma": "W#gamma^{*}",
    "qcd": "QCD Multijet",
    "muon": "muon",
    "electron": "electron"
}

# only for MC samples
variable_dict_mc = {
    "mu": "hnlJet_nominal_isPrompt_MU",
    "gamma": "hnlJet_nominal_isPrompt_PHOTON",
    "tau": "hnlJet_nominal_isPrompt_TAU",
    "e": "hnlJet_nominal_isPrompt_E",
    "uds": "hnlJet_nominal_isUDS",
    "g": "hnlJet_nominal_isG",
    "b": "hnlJet_nominal_isB",
    "c": "hnlJet_nominal_isC",
    "pileup": "hnlJet_nominal_isPU",
    "weight": 'weightNominal',
}

# only for signal samples
variable_dict_hnl = {}

for coupling in [1, 2, 7, 12, 47, 52]:
    variable_dict_hnl[f"weight{coupling}"] = f'weightNominalHNL_{coupling}'

# for all samples (mc bkg/mc sig/data)
variable_dict_all = {
    "tagger": "(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3)*hnlJet_nominal_llpdnnx_ratio_LLP_Q+\
    (nominal_dR_l2j<0.4)*hnlJet_nominal_llpdnnx_ratio_LLP_QMU*subleadingLeptons_isMuon[0]+\
        +(nominal_dR_l2j<0.4)*hnlJet_nominal_llpdnnx_ratio_LLP_QE*subleadingLeptons_isElectron[0]",
    "bdt": "bdt_score_nominal",
    "mass": "nominal_m_llj",
    "met": "nominal_met",
    "mtw": "nominal_mtw",

    "nphoton": "nvetoPhotons",
    "isotropy": "nominal_eventShape_isotropy",
    "circularity": "nominal_eventShape_circularity",
    "aplanarity": "nominal_eventShape_aplanarity",
    "sphericity": "nominal_eventShape_sphericity",

    "boosted": "nominal_dR_l2j<0.4",
    "resolved": "nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3",

    "l2_dxy": "subleadingLeptons_dxy[0]",
    "l2_dxysig": "subleadingLeptons_dxysig[0]",

    "l2_pt": "subleadingLeptons_pt[0]",
    "l2_eta": "subleadingLeptons_eta[0]",
    "l2_iso": "subleadingLeptons_relIso[0]",
}

backgrounds = []
signals = []

process_names = []

ntuple_path = os.path.join(ntuple_path, year)
print(ntuple_path)
cut = weight_dict[category_name]["varexp"]

for process_name, process_dict in samples_dict.items():
    if not (process_name in text_dict.keys()):
        continue
    process = Process(process_name, text_dict[process_name], linecolor='#000000', fillcolor='#000000')

    isMC = True
    if "muon" in process_name or "electron" in process_name:
        isMC=False

    if "HNL" in process_name:
        sample = Sample(process_name, ntuple_path, [process_name+"-"+year], isMC=isMC, year=year, cut=cut)
        process.Add(sample)
    else:
        subprocesses = process_dict[int(year)]
        print(subprocesses)
        for sample_name, sample_list in subprocesses.items():
            print(sample_name)
            sample = Sample(sample_name, ntuple_path, sample_list, isMC=isMC, year=year, oneFile=one_file, cut=cut)
            process.Add(sample)

    vars_to_save = list(variable_dict_all.keys())
    for var, varexp in variable_dict_all.items():
        process.Define(var, varexp)

    if isMC:
        for var, varexp in variable_dict_mc.items():
            process.Define(var, varexp)
        vars_to_save.extend(variable_dict_mc.keys())

    if "HNL" in process_name:
        vars_to_save.extend(variable_dict_hnl.keys())
        for var, varexp in variable_dict_hnl.items():
            process.Define(var, varexp)

    data = process.AsNumpy(*vars_to_save)
    if not os.path.exists(os.path.join(output_dir, year, category_name)):
        os.mkdir(os.path.join(output_dir, year, category_name))

    file_path = os.path.join(output_dir, year, category_name, process_name+".pkl")
    data.to_pickle(file_path)