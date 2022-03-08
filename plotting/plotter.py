import sys
import ROOT
import os
import json
import yaml
import numpy as np
from array import array
import argparse

from histo import Process, Variable, Sample
signal_xsec = 1.

parser = argparse.ArgumentParser()

parser.add_argument("-v", "--variable", action="store", dest="variable")
parser.add_argument("-r", "--region", action="store", dest="region")
parser.add_argument("-y", "--year", action="store", dest="year", default="2016")
parser.add_argument("--region_file", action="store", dest="region_file", default="/home/hep/hsfar/private/limits/histo/config/regions_preselection.json")
parser.add_argument("--variable_file", action="store", dest="variable_file", default="/home/hep/hsfar/private/limits/histo/config/variables.json")
parser.add_argument("--samples_file", action="store", dest="samples_file", default="/home/hep/hsfar/private/limits/histo/config/samples.yml")
parser.add_argument("--ntuple_path", action="store", dest="ntuple_path", default="/vols/cms/hsfar/nanoAOD_friends/13Dec20/")
parser.add_argument("--output_dir", action="store", dest="output_dir", default="plots")
parser.add_argument("--plot_corrections", action="store_true", dest="plotCorrections", default=False)
parser.add_argument("--test", action="store_true", dest="one_file", default=False)

#parser.add_argument("-d", "--data", action="store", dest="data_type")
args = parser.parse_args()
variable_number = int(args.variable)
region_number = int(args.region)
year = str(args.year)
ntuple_path = os.path.join(args.ntuple_path, year)
global weight


with open(os.path.expandvars(args.variable_file)) as json_file:
    lines = json.load(json_file)
    variable_infos = lines[variable_number]
    print(variable_infos)
varname = variable_infos[0].replace('[', '').replace(']', '')

with open(os.path.expandvars(args.region_file)) as json_file:
    lines = json.load(json_file)
    category_infos = lines[region_number]
    category = category_infos[0]
    weight = category_infos[1]
    draw_text = category_infos[2]
    print(category, weight, draw_text)

with open(os.path.expandvars(args.samples_file)) as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)

global data_type
if "el1" in category:
    data_type = "electron"
elif "mu1" in category:
    data_type = "muon"
else:
    print("unknown data")
    sys.exit(-1)
print(year)

def find_xsec(path, xsecs):
    for key, val in xsecs.items():
        if key in path:
            return val

processes = []

text_dict = {
             "wjets": "W+jets",
             "dyjets": "Z/#gamma*+jets",
             "vgamma": "V#gamma*+jets",
             "qcd": "Multijet",
             "topbkg": "t#bar{t}/single-t",
             "muon": "muon",
             "electron": "electron",
             "data": "Data",
             }

if "SR" in category:
    #blinding
    #text_dict["HNL_majorana_all_ctau1p0e02_massHNL4p5_Vall1p016e-03"] = "m_{N} = 4.5 GeV, c#tau_{0} = 10 cm, V_{e}=V_{#mu}"
    text_dict["HNL_majorana_pt20_ctau1p0e-01_massHNL12p0_Vall2p314e-03"] = "m_{N} = 12 GeV, c\\tau_{0} = 0.1 mm, V_{e}=V_{\\mu} (x100)"
    text_dict["HNL_majorana_pt20_ctau1p0e03_massHNL2p0_Vall2p871e-03"] = "m_{N} = 2 GeV, c\\tau_{0} = 1 m, V_{e}=V_{\\mu} (x100)"
    text_dict.pop("data")
    text_dict.pop("muon")
    text_dict.pop("electron")


for process_name in text_dict.keys():
    if process_name not in samples_dict.keys():
        continue
    process = samples_dict[process_name]
    isMC=True
    if "muon" in process_name or "electron" in process_name:
        if process_name!=data_type:
            continue
        else:
            process_name = "data"
            isMC = False

    if "HNL" in process_name:
        process = Process("HNL", text_dict[process_name], linecolor='#bd0000')
        process.Add(Sample(process_name, ntuple_path, ["{}-{}".format(process_name, year)], year=year, limits=False, cut=weight))

    else:
        subprocesses = process[int(year)]
        process = Process(process_name, text_dict[process_name])
        for sample_name, sample_list in subprocesses.items():
            print(sample_name)
            sample = Sample(sample_name, ntuple_path, sample_list, isMC=isMC, year=year, oneFile=args.one_file, cut=weight)
            process.Add(sample)
        
    process.Define(varname.replace("-", "m").replace("(","_").replace(")","_").replace("/", "over"), variable_infos[0])
    processes.append(process)

variable_infos[0] = varname
variable = Variable(*variable_infos, corrections=args.plotCorrections)
for process in processes:
    print(process.name)
    if "HNL" in process.name:
        hist = process.Histo1D(variable.args, varname.replace("-", "m").replace("(","_").replace(")","_").replace("/", "over"), weight="weightNominalHNL_7")
        hist.Scale(100.)
        variable.Add(hist, process.title, isSignal=True, isData=False)
    else:
        if process.name == "data":
            isData = True
        else:
            isData = False
        hist = process.Histo1D(variable.args, varname.replace("-", "m").replace("(","_").replace(")","_").replace("/", "over"))
        variable.Add(hist, process.title, isSignal=False, isData=isData)

        if args.plotCorrections and not isData:
            histUp = process.Histo1D(variable.args, varname.replace("-", "m").replace("(","_").replace(")","_").replace("/", "over"), weight="weightNominalCorrectedUp")
            variable.Add(histUp, process.title, isSignal=False, isData=isData, correction="up")

variable.Draw(category, "hist", draw_text, year=year, output_dir=args.output_dir)
