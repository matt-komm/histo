import ROOT
import os
import json

procs = ["HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03-2016", "WToLNu_0J_13TeV-amcatnloFXFX-pythia8-ext1-2016", "WToLNu_1J_13TeV-amcatnloFXFX-pythia8-2016", "WToLNu_2J_13TeV-amcatnloFXFX-pythia8-ext4-2016", "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016" , "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016"]
ntuple_path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21/2016"


for proc in procs:
    print(proc)
    file_list = ROOT.std.vector('string')()

    for f in os.listdir(os.path.join(ntuple_path, proc)):
        file_list.push_back(os.path.join(ntuple_path, proc, f))

    rdf = ROOT.RDataFrame("Friends", file_list)

    rdf = rdf.Filter("nselectedJets_nominal>0 and nselectedJets_nominal<5", "njet")
    rdf = rdf.Filter("MET_filter and nominal_met<100", "met")
    rdf = rdf.Filter("nominal_dR_l2j<1.3", "dR_l2j")
    rdf = rdf.Filter("dilepton_mass<80 and dilepton_mass>20", "mll")
    rdf = rdf.Filter("nominal_m_llj<90 and nominal_m_llj>70", "mllj")

    rdf.Report().Print()
