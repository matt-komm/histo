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

categories_file = os.path.expandvars("$HISTO_BASE_ENV/config/categories_2l_inclusive.json")
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

backgrounds = []
signals = []

prediction_dict = {}
observed_dict = {}
prediction_error_dict = {}
observed_error_dict = {}
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


    tagger_score = "(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3)*hnlJet_nominal_llpdnnx_ratio_LLP_Q+\
    (nominal_dR_l2j<0.4)*hnlJet_nominal_llpdnnx_ratio_LLP_QMU*subleadingLeptons_isMuon[0]+\
        +(nominal_dR_l2j<0.4)*hnlJet_nominal_llpdnnx_ratio_LLP_QE*subleadingLeptons_isElectron[0]"

    mass_variable = "nominal_m_llj"

    if isMC:
        process.Define("mu", "hnlJet_nominal_isPrompt_MU")
        process.Define("gamma", "hnlJet_nominal_isPrompt_PHOTON")
        process.Define("tau", "hnlJet_nominal_isPrompt_TAU")
        process.Define("e", "hnlJet_nominal_isPrompt_E")
        process.Define("uds", "hnlJet_nominal_isUDS")
        process.Define("g", "hnlJet_nominal_isG")
        process.Define("b", "hnlJet_nominal_isB")
        process.Define("c", "hnlJet_nominal_isC")
        process.Define("pileup", "hnlJet_nominal_isPU")

    process.Define("ntracks", "hnlJet_nominal_ncpf")
    process.Define("tagger_score", tagger_score)
    process.Define("bdt_score_nominal", "bdt_score_nominal")
    process.Define("mass", mass_variable)
    process.Define("met", "nominal_met")
    process.Define("mtw", "nominal_mtw")

    process.Define("nphoton", "nvetoPhotons")
    process.Define("isotropy", "nominal_eventShape_isotropy")
    process.Define("circularity", "nominal_eventShape_circularity")
    process.Define("aplanarity", "nominal_eventShape_aplanarity")
    process.Define("sphericity", "nominal_eventShape_sphericity")

    process.Define("weight", 'weightNominal')
    if "HNL" in process.name:
        for coupling in [1, 2, 7, 12, 47, 52]:
            process.Define(f"weight{coupling}", f'weightNominalHNL_{coupling}')

    process.Define("boosted", "nominal_dR_l2j<0.4")
    process.Define("resolved", "nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3")

    process.Define("dxy", "subleadingLeptons_dxy[0]")
    process.Define("dxysig", "subleadingLeptons_dxysig[0]")

    process.Define("l2_pt", "subleadingLeptons_pt[0]")
    process.Define("l2_eta", "subleadingLeptons_eta[0]")
    process.Define("l2_iso", "subleadingLeptons_relIso[0]")


    vars_to_save = ["mass", "tagger_score", "boosted", "resolved", "bdt_score_nominal", "dxy", "dxysig", "weight", "dilepton_mass", "dilepton_dPhimin", "ntracks", "nphoton", "hnlJet_nominal_pt", "hnlJet_nominal_eta", "l2_pt", "l2_eta", "l2_iso", "isotropy", "circularity", "aplanarity", "sphericity", "met", "mtw"]
    jet_flavours = ["mu", "gamma", "tau", "e", "uds", "g", "b", "c", "pileup"]
    if isMC:
        vars_to_save.extend(jet_flavours)

    if "HNL" in process_name:
        vars_to_save.extend(["weight1", "weight2", "weight7", "weight12", "weight47", "weight52"])
    data = process.AsNumpy(*vars_to_save)
    #data['merged'] = data['merged'].astype(bool)
    #data['resolved'] = data['resolved'].astype(bool)
    #data['tagger_score'] = data['tagger_score'].astype(float)
    bkg_yield = np.sum(data['weight'])
    print(f"yield {np.sum(data['weight'])}")
    print(os.path.join(output_dir, year, category_name))
    if not os.path.exists(os.path.join(output_dir, year, category_name)):
        os.mkdir(os.path.join(output_dir, year, category_name))

    file_path = os.path.join(output_dir, year, category_name, process_name+".pkl")
    data.to_pickle(file_path)

    # if isMC:
    #     flavour_sum = 0
    #     for jet_flavour in jet_flavours:

    #         flavour_yield = np.sum(data.query(f"{jet_flavour}==1")['weight'])

    #         print(f"fraction of {jet_flavour}: {flavour_yield/bkg_yield}")
    #         flavour_sum += flavour_yield/bkg_yield
    #     print(f"Total sum of all flavours: {flavour_sum}")
