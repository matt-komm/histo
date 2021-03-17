import ROOT
import os
import yaml
import json
import argparse
import classes.process
import classes.variable
import classes.sample


# make a new hist for categorised events passing the cuts and write it to file
def write_hist(hist_nano, category_dict, name, isMC=True):
    index_new = 0
    hist_limits = ROOT.TH1D("hist", "hist", len(category_dict), 0, len(category_dict))

    for index, category in category_dict.items():
        index_new += 1
        hist_content = hist_nano.GetBinContent(hist_nano.FindBin(index))
        hist_error = hist_nano.GetBinError(hist_nano.FindBin(index))

        print(index, index_new, category, hist_content)
        hist_limits.SetBinContent(index_new, hist_content)
        hist_limits.SetBinError(index_new, hist_error)
        hist_limits.GetXaxis().SetBinLabel(index_new, category)

    if isMC:
        removeNegEntries(hist_limits)
    hist_limits.SetTitle(name)
    hist_limits.SetName(name)
    hist_limits.SetDirectory(root_file)
    hist_limits.Write()


def removeNegEntries(hist):
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


# make histograms per year, process
parser = argparse.ArgumentParser()
parser.add_argument("--year",default="2016")
parser.add_argument("--proc", default="qcd")
parser.add_argument("--category", default="mumu_OS")
parser.add_argument("--region", default="D")
parser.add_argument("--ntuple_path", default="/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/22Dec20/") #03Dec20
parser.add_argument("--data", action="store_true", default=False)
parser.add_argument("--one_file", action="store_true", dest="oneFile", default=False)

args = parser.parse_args()
year = args.year
proc = args.proc
macroCategory_name = args.category
ntuple_path = os.path.join(args.ntuple_path, year)
region = args.region
oneFile = args.oneFile
isData = args.data
isMC = not isData

lumi = {"2016": 35.92, "2017": 41.53, "2018": 59.68}

with open("samples.yaml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)
    subprocesses = samples_dict[proc]

#####################################
### Various configurations go here

# Systematic uncertainties
systDict = {}
systDict["IsoMuTrigger_weight_trigger"] = "trigger"
systDict["tightMuon_weight_iso"] = "tight_muon_iso"
systDict["tightMuon_weight_id"] = "tight_muon_id"
systDict["tightElectron_weight_id"] = "tight_electron_id"
systDict["tightElectron_weight_reco"] = "tight_electron_reco"
systDict["looseElectrons_weight_reco"] = "loose_electron_reco"
systDict["puweight"] = "pu"

systematics_shapes = ["nominal", "jesTotalUp", "jesTotalDown", "jerUp", "jerDown"]

####################################

# For data-driven background estimation
def abcd_cut(region="D", syst="nominal", isData=False, single=False):
    cut = ""
    if single:
        if region == "A" or region == "B":
            cut += f'((category_simplified_{syst}_single_index==2 and (category_simplified_{syst}_llpdnnx_m_llj>90. or category_simplified_{syst}_llpdnnx_m_llj<60.)) \
                or (category_simplified_{syst}_single_index==1 and (category_simplified_{syst}_llpdnnx_m_lljj>95. or category_simplified_{syst}_llpdnnx_m_lljj<75.)))'
        elif region == "C" or region == "D":
            cut += f'((category_simplified_{syst}_single_index==2 and category_simplified_{syst}_llpdnnx_m_llj<90. and category_simplified_{syst}_llpdnnx_m_llj>70.) \
                or (category_simplified_{syst}_single_index==1 and category_simplified_{syst}_llpdnnx_m_lljj<95. and category_simplified_{syst}_llpdnnx_m_lljj>75.))'
    else: #dilep
        if region == "A" or region == "B":
            cut += f'(category_simplified_{syst}_llpdnnx_m_llj>90. or category_simplified_{syst}_llpdnnx_m_llj<70.)'
        elif region == "C" or region == "D":
            cut += f'(category_simplified_{syst}_llpdnnx_m_llj<90.) and (category_simplified_{syst}_llpdnnx_m_llj>70.)'


    if single:
        if not isData:
            if region == "A" or region == "D":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max2nd>30 and category_simplified_{syst}_single_index==1)'
                cut += f' or (category_simplified_{syst}_llpdnnx_max>500 and category_simplified_{syst}_single_index==2))'
            elif region == "B" or region == "C":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max2nd<30 and category_simplified_{syst}_single_index==1)'
                cut += f' or (category_simplified_{syst}_llpdnnx_max>30 and category_simplified_{syst}_llpdnnx_max<500 and category_simplified_{syst}_single_index==2))'
        else:
            if region == "A" or region == "D":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max2nd>15 and category_simplified_{syst}_llpdnnx_max2nd<20 and category_simplified_{syst}_single_index==1)'
                cut += f' or (category_simplified_{syst}_llpdnnx_max>100 and category_simplified_{syst}_llpdnnx_max<300 and category_simplified_{syst}_single_index==2))'
            elif region == "B" or region == "C":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max2nd<15 and category_simplified_{syst}_single_index==1)'
                cut += f' or (category_simplified_{syst}_llpdnnx_max>30 and category_simplified_{syst}_llpdnnx_max<100 and category_simplified_{syst}_single_index==2))'

    else: #dilep
        if not isData:
            if region == "A" or region == "D":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max>100 and category_simplified_{syst}_index==1) \
                          or (category_simplified_{syst}_llpdnnx_max>300 and category_simplified_{syst}_index==2))'
            elif region == "B" or region == "C":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max>1 and category_simplified_{syst}_llpdnnx_max<100 and category_simplified_{syst}_index==1) \
                          or (category_simplified_{syst}_llpdnnx_max>1 and category_simplified_{syst}_llpdnnx_max<300 and category_simplified_{syst}_index==2))'
        else:
            if region == "A" or region == "D":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max>10 and category_simplified_{syst}_llpdnnx_max<30 and category_simplified_{syst}_index==1) \
                          or (category_simplified_{syst}_llpdnnx_max>30 and category_simplified_{syst}_llpdnnx_max<100 and category_simplified_{syst}_index==2))'
            elif region == "B" or region == "C":
                cut += f' and ((category_simplified_{syst}_llpdnnx_max>1 and category_simplified_{syst}_llpdnnx_max<10 and category_simplified_{syst}_index==1) \
                          or (category_simplified_{syst}_llpdnnx_max>1 and category_simplified_{syst}_llpdnnx_max<30 and category_simplified_{syst}_index==2))'

    return cut

#####################################

# couplings to consider
couplings = [2, 7, 12, 47, 52, 67]

with open('categories.json', 'r') as fp:
    macroCategory_dict = json.load(fp)

macroCategory_cut = macroCategory_dict[macroCategory_name]

diLeptonCategory_dict = {}
diLeptonCategory_dict[1] = "q"
diLeptonCategory_dict[2] = "ql"

singleLeptonCategory_dict = {}
singleLeptonCategory_dict[1] = "qq"
singleLeptonCategory_dict[2] = "q"
#####################################


if "HNL" in proc:
    process = classes.process.Process("HNL", proc)
    process.add(classes.sample.Sample(proc, ntuple_path, ["{}-{}".format(proc, year)], year=year, limits=True))
else:
    process = classes.process.Process(proc, proc)
    subprocesses = subprocesses[int(year)]
    for sample_name, sample_list in subprocesses.items():
        print(sample_name)
        sample = classes.sample.Sample(sample_name, ntuple_path, sample_list, year=year, oneFile=oneFile, isMC=isMC)
        process.add(sample)


if isMC:
    for systName, abrv in systDict.items():
        for variation in ["Up", "Down"]:
            if "HNL" in process.name:
                for coupling in range(2, 68):
                    process.Define("weightHNL_{}_{}{}".format(coupling, abrv, variation), "weightNominalHNL_{}/{}*{}".format(coupling, systName+"_nominal", systName+"_"+variation.lower()))
            else:
                process.Define("weight_{}{}".format(abrv, variation), "weightNominal/{}*{}".format(systName+"_nominal", systName+"_"+variation.lower()))

category_variable_nominal = "category_simplified_nominal_index"  

is_single_lep = True if "single" in macroCategory_name else False
if is_single_lep:
    # Define custom variable to perform single-lep categorisation into 1 and 2-tag
    for systName in systematics_shapes:
        if isData:
            if systName != "nominal":
                continue
        process.Define(f"category_simplified_{systName}_single_index", f"(category_simplified_{systName}_llpdnnx_max2nd<10)*2+(category_simplified_{systName}_llpdnnx_max2nd>10)*1")
    category_variable_nominal = "category_simplified_nominal_single_index"

# create root file with nominal value histogram and various systematic variations
# to be used with Combine Harvester
root_file = ROOT.TFile.Open("hists/{}_{}_{}_{}.root".format(proc, args.category, region, year), "RECREATE")
print("The category name and cut are:", macroCategory_name, macroCategory_cut)
root_file.cd()
root_file.mkdir(macroCategory_name+"_"+region)
root_file.cd(macroCategory_name+"_"+region)

if is_single_lep:
    category_dict = singleLeptonCategory_dict
else:
    category_dict = diLeptonCategory_dict


coupling = 1
while coupling < 67:
    # different scenarios
    # Need to calculate yield per coupling
    if "HNL" in process.name:
        coupling += 1
        if coupling not in couplings:
            continue
    else:
        coupling = 68


    # variations with changing shape and constant weight
    if isMC:
        for systName in systematics_shapes:
            if "HNL" in process.name:
                name = "{}_coupling_{}".format(process.name, coupling)
            else:
                name = process.name

            # Need to be changed for hnl
            if "HNL" in process.name:
                weight = "weightNominalHNL_{}".format(coupling)
            else:
                weight = "weightNominal"
            
            # add name for variations
            if systName != "nominal":
                name += "_"+systName

            # Systematic variation -- replace nominal by systematic in all cuts
            cut = macroCategory_cut.replace("nominal", systName)
            cut += " and {}".format(abcd_cut(syst=systName, region=region, single=is_single_lep))

            # Variable used to categorise the events. replace nominal by systematic
            category_variable = category_variable_nominal.replace("nominal", systName)

            # read in hist from nanoAOD friends
            print(systName, category_variable, cut, weight)

            hist_nano = process.Histo1D((category_variable, category_variable, 7, -2, 5), category_variable, cut=cut, weight=weight)
            hist_nano = hist_nano.Clone()
            write_hist(hist_nano, category_dict, name)


        # variations with constant shape but changing weight
        for systName, abrv in systDict.items():
            for variation in ["Up", "Down"]:
                if "HNL" in process.name:
                    name = "{}_coupling_{}_{}{}".format(process.name, coupling, abrv, variation)
                else:
                    name = "{}_{}{}".format(process.name, abrv, variation)

                if "HNL" in process.name:
                    weight = "weightHNL_{}_{}{}".format(coupling, abrv, variation)
                else:
                    weight = "weight_{}{}".format(abrv, variation)

                # only nominal cut will be applied
                cut = macroCategory_cut
                cut += " and {}".format(abcd_cut(syst="nominal", region=region, single=is_single_lep))
                category_variable = category_variable_nominal
                print(systName, category_variable, cut, weight)

                # read in hist from nanoAOD friends
                hist_nano = process.Histo1D((category_variable, category_variable, 7, -2, 5), category_variable, cut=cut, weight=weight)
                hist_nano = hist_nano.Clone()
                write_hist(hist_nano, category_dict, name)
    else:
        name = process.name
        cut = macroCategory_cut
        cut += " and {}".format(abcd_cut(syst="nominal", region=region, isData=isData))
        category_variable = category_variable_nominal
        print(category_variable, cut)

        hist_nano = process.Histo1D((category_variable, category_variable, 7, -2, 5), category_variable, cut=cut, weight="weightNominal")
        hist_nano = hist_nano.Clone()
        write_hist(hist_nano, category_dict, "data", isMC=False)

root_file.Close()
