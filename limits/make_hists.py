import ROOT
import os
import yaml
import json
import argparse
from histo import Process, Sample

def write_hist(hist_nano: ROOT.TH1D, category_dict: dict, name: str, isMC: bool=True):
    """
    Make a new histogram for categorised events passing the cuts and write it to file

    Args:
        hist_nano     : nanoAOD histogram
        category_dict : dictionary of category cuts
        name : category name
        isMC : flag for MC/data
    """

    index_new = 0
    hist_limits = ROOT.TH1D("hist", "hist", len(category_dict), 0, len(category_dict))
    print(f"Writing hist {name}")
    for index, category in category_dict.items():
        index_new += 1
        hist_content = hist_nano.GetBinContent(hist_nano.FindBin(index))
        hist_error = hist_nano.GetBinError(hist_nano.FindBin(index))
        print(index, index_new, category, hist_content)
        hist_limits.SetBinContent(index_new, hist_content)
        hist_limits.SetBinError(index_new, hist_error)
        hist_limits.GetXaxis().SetBinLabel(index_new, category)

    if isMC:
        remove_neg_entries(hist_limits)
    hist_limits.SetTitle(name)
    hist_limits.SetName(name)
    hist_limits.SetDirectory(root_file)
    hist_limits.Write()


def remove_neg_entries(hist: ROOT.TH1D):
    """ Removes negative and/or small entries from a histogram """

    alpha = 1. - 0.6827
    upErr = ROOT.Math.gamma_quantile_c(alpha/2,1,1)
    avgWeight = hist.Integral()/hist.GetEntries() if hist.GetEntries()>0 else -1
    #print "weight",avgWeight
    for ibin in range(hist.GetNbinsX()):
        c = hist.GetBinContent(ibin+1)
        if c<10**-4:
            hist.SetBinContent(ibin+1,10**-3)
            #note: in case of 0 entries the uncertainty is also small
            #(this is not the case with negative events)
            if hist.GetBinError(ibin+1)<10**-4 and avgWeight>0:
                #set uncertainties for empy bins
                #https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
                hist.SetBinError(ibin+1,upErr*avgWeight)
            else:
                hist.SetBinError(ibin+1,10**-4)
        #print "bin%2i, %.1f+-%.1f (+-%.1f%%)"%(ibin,c,hist.GetBinError(ibin+1),100.*hist.GetBinError(ibin+1)/c if c>0 else -1)


def mass_cut(delta_m:float=5., region:str="D", syst:str="nominal", single_lepton=False) -> str:
    """ 
    Returns mass window cut repending on region

    Args:
    delta_m : mW+-delta_m used to define signal region band
    region: A, B, C, D
    syst: systematic variation
    single_lepton: deprecated

    Returns:
    mass cut str 
    """
    mW=80.
    m_upper = mW + delta_m
    m_lower = mW - delta_m

    if region not in ["A", "B", "C", "D"]:
        raise ValueError(f"Invalid region {region} selected")

    if region == "A" or region == "C":
        return  f'({syst}_m_llj<{m_lower}) '
        #return  f'({syst}_m_llj>{m_upper} or {syst}_m_llj<{m_lower}) '
    elif region == "B" or region == "D":
        return  f'({syst}_m_llj<{m_upper} and {syst}_m_llj>{m_lower}) '
    else:
        return ""

def bdt_cut(bdt_threshold: float, syst:str="nominal") -> str:
     return f'(bdt_score_{syst} > {bdt_threshold})'

def tagger_cut(tagger_threshold: float, lower_threshold: float=0.1, region:str="D", syst:str="nominal") -> str:
    """ Returns tagger score cut repending on region 
    Args:
    tagger_threshold : define signal region
    lower_threshold: lower tagger score threshold
    region: A, B, C, D
    syst: systematic variation

    Returns:
    tagger cut str 
    """
    if region not in ["A", "B", "C", "D"]:
        raise ValueError(f"Invalid region {region} selected")

    if tagger_threshold < lower_threshold:
        raise ValueError("Inconsistent tagger thresholds")

    if region == "A" or region == "B":
        return f'(tagger_score_{syst} > {lower_threshold} and tagger_score_{syst} < {tagger_threshold})'
    elif region == "C" or region == "D":
        return f'(tagger_score_{syst} > {tagger_threshold})'


def tagger_compound_variable(syst:str="nominal", single_lepton=False) -> str:
    """ Compound tagger variable to facilitate categorisation """
    # if single_lepton:
    #     return f"hnlJet_{syst}_llpdnnx_ratio_LLP_Q"
    # else:
    return f"({syst}_dR_l2j>0.4 and {syst}_dR_l2j<1.3)*hnlJet_{syst}_llpdnnx_ratio_LLP_Q+\
        ({syst}_dR_l2j<0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_QMU*subleadingLeptons_isMuon[0]+\
        +({syst}_dR_l2j<0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_QE*subleadingLeptons_isElectron[0]"

def make_hists(process, systematics_shapes, systematics_rates, cut_nominal, category_variable_nominal, thresholds, region, coupling=None):

    hists = {}

    def make_hist(process, category_variable, thresholds, weight, cut, region, syst="nominal"):
        thresholds_merged_prompt = thresholds["prompt"]["merged"]
        thresholds_merged_medium = thresholds["medium"]["merged"]
        thresholds_merged_displaced = thresholds["displaced"]["merged"]
        thresholds_resolved_prompt = thresholds["prompt"]["resolved"]
        thresholds_resolved_medium = thresholds["medium"]["resolved"]
        thresholds_resolved_displaced = thresholds["displaced"]["resolved"]

        #print(f"Tagger thresholds, merged: {threshold_merged}, resolved: {threshold_resolved}")
        #print(f"DeltaM, merged: {deltam_merged}, resolved: {deltam_resolved}")
        mass_cut_merged_prompt = mass_cut(delta_m=thresholds_merged_prompt[1], region=region, syst=syst)
        mass_cut_merged_medium = mass_cut(delta_m=thresholds_merged_medium[1], region=region, syst=syst)
        mass_cut_merged_displaced = mass_cut(delta_m=thresholds_merged_displaced[1], region=region, syst=syst)

        mass_cut_resolved_prompt = mass_cut(delta_m=thresholds_resolved_prompt[1], region=region, syst=syst)
        mass_cut_resolved_medium = mass_cut(delta_m=thresholds_resolved_medium[1], region=region, syst=syst)
        mass_cut_resolved_displaced = mass_cut(delta_m=thresholds_resolved_displaced[1], region=region, syst=syst)

        tagger_cut_merged_prompt = tagger_cut(thresholds_merged_prompt[0], region=region, syst=syst)
        tagger_cut_merged_medium = tagger_cut(thresholds_merged_medium[0], region=region, syst=syst)
        tagger_cut_merged_displaced = tagger_cut(thresholds_merged_displaced[0], region=region, syst=syst)

        tagger_cut_resolved_prompt = tagger_cut(thresholds_resolved_prompt[0], region=region, syst=syst)
        tagger_cut_resolved_medium = tagger_cut(thresholds_resolved_medium[0], region=region, syst=syst)
        tagger_cut_resolved_displaced = tagger_cut(thresholds_resolved_displaced[0], region=region, syst=syst)

        bdt_cut_merged_prompt = bdt_cut(thresholds_merged_prompt[2], syst=syst)
        bdt_cut_merged_medium = bdt_cut(thresholds_merged_medium[2], syst=syst)
        bdt_cut_merged_displaced = bdt_cut(thresholds_merged_displaced[2], syst=syst)

        bdt_cut_resolved_prompt = bdt_cut(thresholds_resolved_prompt[2], syst=syst)
        bdt_cut_resolved_medium = bdt_cut(thresholds_resolved_medium[2], syst=syst)
        bdt_cut_resolved_displaced = bdt_cut(thresholds_resolved_displaced[2], syst=syst)

        cut += f" and ("
        cut += f"({category_variable}==1 and {mass_cut_merged_prompt} and {tagger_cut_merged_prompt} and {bdt_cut_merged_prompt}) or "
        cut += f"({category_variable}==2 and {mass_cut_merged_medium} and {tagger_cut_merged_medium} and {bdt_cut_merged_medium}) or "
        cut += f"({category_variable}==3 and {mass_cut_merged_displaced} and {tagger_cut_merged_displaced} and {bdt_cut_merged_displaced}) or "
        cut += f"({category_variable}==4 and {mass_cut_resolved_prompt} and {tagger_cut_resolved_prompt} and {bdt_cut_resolved_prompt}) or "
        cut += f"({category_variable}==5 and {mass_cut_resolved_medium} and {tagger_cut_resolved_medium} and {bdt_cut_resolved_medium}) or "
        cut += f"({category_variable}==6 and {mass_cut_resolved_displaced} and {tagger_cut_resolved_displaced} and {bdt_cut_resolved_displaced})"
        cut += f")"

        print(syst, category_variable, cut, weight)

        hist_nano = process.Histo1D((category_variable, category_variable, 6, 0.5, 6.5), category_variable, cut=cut, weight=weight)
        hist_nano = hist_nano.Clone()
        return hist_nano
    if systematics_rates is not None:
    # variations with constant shape but changing weight
        for syst, abrv in systematics_rates.items():
            for variation in ["Up", "Down"]:
                if "HNL" in process.name:
                    name = f"{process.name}_coupling_{coupling}_{abrv}{year}{variation}"
                    weight = f"weightHNL_{coupling}_{abrv}{variation}"
                else:
                    name = f"{process.name}_{abrv}{year}{variation}"
                    weight = f"weight_{abrv}{variation}"
                hists[name] = make_hist(process, category_variable_nominal, thresholds, weight, cut_nominal, region, syst="nominal")

    if systematics_shapes is not None:
        for syst in systematics_shapes:
            # add name for variations
            variations = ["Up", "Down"] if "nominal" not in syst else [""]
            syst_hack = syst.replace("nominal", "")
            for variation in variations:
                if "HNL" in process.name:
                    if "nominal" in syst:
                        name = f"{process.name}_coupling_{coupling}"
                    else:
                        name = f"{process.name}_coupling_{coupling}_{syst_hack}{year}{variation}"

                    weight = f"weightNominalHNL_{coupling}"
                else:
                    if "nominal" in syst:
                        name = f"{process.name}"
                    else:
                        name = f"{process.name}_{syst_hack}{year}{variation}"
                    weight = "weightNominal"

                syst_var_name = f"{syst}{variation}"

                # Systematic variation -- replace nominal by systematic in all cuts
                cut = cut_nominal.replace("nominal", syst_var_name)
                cut = cut.replace("nselectedJets_unclEnUp", "nselectedJets_nominal") #Hack!
                cut = cut.replace("nselectedJets_unclEnDown", "nselectedJets_nominal") #Hack!
                category_variable = category_variable_nominal.replace("nominal", syst_var_name)
                # read in hist from nanoAOD friends
                hists[name] = make_hist(process, category_variable, thresholds, weight, cut, region, syst=syst_var_name)

    return hists

# make histograms per year, process
parser = argparse.ArgumentParser()
parser.add_argument("--year",default="2016")
parser.add_argument("--proc", default="wjets")
parser.add_argument("--category", default="mumu_OS")
parser.add_argument("--region", default="D")
parser.add_argument("--ntuple_path", default="/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/28May21")
parser.add_argument("--output_path", default="hists")
parser.add_argument("--data", action="store_true", default=False)
parser.add_argument("--test", action="store_true", dest="oneFile", default=False)

args = parser.parse_args()
print(vars(args))

year = args.year
proc = args.proc
category_name = args.category
ntuple_path = os.path.join(f"{args.ntuple_path}", year)
region = args.region
oneFile = args.oneFile
isData = args.data
isMC = not isData
output_path = args.output_path

with open("../config/samples.yml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)
    subprocesses = samples_dict[proc]

#####################################
### Various configurations go here

# Systematic uncertainties
systematics_rates = {}
systematics_rates["IsoMuTrigger_weight_trigger"] = "trigger"
systematics_rates["tightMuons_weight_iso"] = "tight_muon_iso"
systematics_rates["tightMuons_weight_id"] = "tight_muon_id"
systematics_rates["tightElectrons_weight_id"] = "tight_electron_id"
systematics_rates["tightElectrons_weight_reco"] = "tight_electron_reco"
systematics_rates["looseElectrons_weight_reco"] = "loose_electron_reco"
systematics_rates["puweight"] = "pu"
systematics_rates["hnlJet_track_weight_adapted"] = "displaced_track"

systematics_shapes = ["nominal", "jesTotal", "jer", "unclEn"]
####################################

# couplings to consider
#couplings = range(2, 68)
couplings = [1, 2, 7, 12, 47, 52]
#couplings = [7]

category_file = '../config/categories_2l_inclusive.json'
threshold_file = f'../config/coordsBestThresholdsIncBdt_{year}.json'

with open(category_file, 'r') as fp:
    categories_2l = json.load(fp)

with open(threshold_file, 'r') as fp:
    threshold_dict = json.load(fp)

category_cut = categories_2l[category_name]["varexp"]
thresholds = {}
thresholds["prompt"] = threshold_dict[category_name+"_prompt"]
thresholds["medium"]= threshold_dict[category_name+"_medium"]
thresholds["displaced"]= threshold_dict[category_name+"_displaced"]


dilepton_category_dict = {}
dilepton_category_dict[1] = "ql, dxysig<1"
dilepton_category_dict[2] = "ql, 1<dxysig<10"
dilepton_category_dict[3] = "ql, dxysig>10"
dilepton_category_dict[4] = "q, dxysig<1"
dilepton_category_dict[5] = "q, 1<dxysig<10"
dilepton_category_dict[6] = "q, dxysig>10"

# Process configuration
if "HNL" in proc:
    process = Process("HNL", proc)
    process.Add(Sample(proc, ntuple_path, ["{}-{}".format(proc, year)], year=year, limits=True))
else:
    process = Process(proc, proc)
    subprocesses = subprocesses[int(year)]
    for sample_name, sample_list in subprocesses.items():
        print(sample_name)
        sample = Sample(sample_name, ntuple_path, sample_list, year=year, oneFile=oneFile, isMC=isMC)
        process.Add(sample)

# Event weights: MC only
if isMC:
    for syst, abrv in systematics_rates.items():
        for variation in ["Up", "Down"]:
            if "HNL" in process.name:
                for coupling in couplings:
                    process.Define("weightHNL_{}_{}{}".format(coupling, abrv, variation), "weightNominalHNL_{}/{}*{}".format(coupling, syst+"_nominal", syst+"_"+variation.lower()))
            else:
                process.Define("weight_{}{}".format(abrv, variation), "weightNominal/{}*{}".format(syst+"_nominal", syst+"_"+variation.lower()))
    for syst in systematics_shapes:
        for variation in ["Up", "Down"]:
            if "nominal" in syst:
                variation = ""
            # Define resolved and merged categories & tagger score variable & mass cuts
            #process.Define(f"category_{syst}_index", f"1.*({syst}_dR_l2j<0.4) \
                                                    #+ 2.*({syst}_dR_l2j>0.4 and {syst}_dR_l2j<1.3)")
            process.Define(f"category_{syst+variation}_index", f"1.*({syst+variation}_dR_l2j<0.4 and subleadingLeptons_dxysig[0]<1.) \
                                                    + 2.*({syst+variation}_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>1. and subleadingLeptons_dxysig[0]<10.)\
                                                    + 3.*({syst+variation}_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>10.) \
                                                    + 4.*({syst+variation}_dR_l2j>0.4 and {syst+variation}_dR_l2j<1.3 and subleadingLeptons_dxysig[0]<1.) \
                                                    + 5.*({syst+variation}_dR_l2j>0.4 and {syst+variation}_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>1. and subleadingLeptons_dxysig[0]<10.) \
                                                    + 6.*({syst+variation}_dR_l2j>0.4 and {syst+variation}_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>10.)")

            process.Define(f"tagger_score_{syst+variation}", tagger_compound_variable(syst+variation, single_lepton=False))

else:
    process.Define(f"category_nominal_index", f"1.*(nominal_dR_l2j<0.4 and subleadingLeptons_dxysig[0]<1.) \
                                                + 2.*(nominal_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>1. and subleadingLeptons_dxysig[0]<10.)\
                                                + 3.*(nominal_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>10.) \
                                                + 4.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3 and subleadingLeptons_dxysig[0]<1.) \
                                                + 5.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>1. and subleadingLeptons_dxysig[0]<10.) \
                                                + 6.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>10.)")
    process.Define(f"tagger_score_nominal", tagger_compound_variable(syst="nominal", single_lepton=False))


# create root file with nominal value histogram and various systematic variations
# to be used with Combine Harvester
root_file = ROOT.TFile.Open(os.path.join(output_path, f"{proc}_{args.category}_{region}_{year}.root"), "RECREATE")
print("The category name and cut are:", category_name, category_cut)
root_file.cd()
root_file.mkdir(category_name+"_"+region)
root_file.cd(category_name+"_"+region)

category_dict = dilepton_category_dict
category_variable_nominal = "category_nominal_index"

coupling = 1
while coupling < 67:
    # different scenarios
    # Need to calculate yield per coupling
    if "HNL" in process.name:
        if coupling not in couplings:
            continue
    else:
        coupling = 68

    if isMC:
        hists = make_hists(process, systematics_shapes, systematics_rates, category_cut, category_variable_nominal, thresholds, region, coupling=coupling)
        for name, hist in hists.items():
            write_hist(hist, category_dict, name, isMC=True)
    else:
        hists = make_hists(process, ["nominal"], None, category_cut, category_variable_nominal, thresholds, region, coupling=coupling)
        for name, hist in hists.items():
            write_hist(hist, category_dict, "data", isMC=False)
    coupling += 1


root_file.Close()
