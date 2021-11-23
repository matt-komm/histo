import ROOT
import numpy as np
import yaml
import argparse
import os
import random
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument("--year",type=int,default=2016, choices=[2016,2017,2018])
parser.add_argument("--plot", default="mllj_SR")
parser.add_argument("-o","--output", dest='output', default="hists")
parser.add_argument("--syst", default="nominal")
args = parser.parse_args()

#/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21/2016

yields = json.load(open("/vols/cms/LLP/yields_201117/"+str(args.year)+"/eventyields.json"))
yieldsHNL = json.load(open("/vols/cms/LLP/yields_201117/"+str(args.year)+"/eventyieldsHNL.json"))

xsecs = json.load(open("/vols/cms/LLP/xsec.json"))
hnlXsecs =  json.load(open("/vols/cms/LLP/gridpackLookupTable.json"))
hnlFilter = json.load(open("/vols/cms/LLP/filterTable.json"))

lumiPerYear = {2016: 35.92, 2017: 41.53, 2018: 59.68}

#print (hnlXsecs)




    

def njetsVar(syst="nominal"):
    if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown']:
        return "nselectedJets_"+syst
    else:
        return "nselectedJets_nominal"

def nfwdjetsVar(syst="nominal"):
    if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown']:
        return "nselectedFwdJets_"+syst
    else:
        return "nselectedFwdJets_nominal"

def metVar(syst="nominal"):
    if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown', 'unclEnUp', 'unclEnDown']:
        return syst+"_met"
    else:
        return "nominal_met"

def bdtVar(syst="nominal"):
    if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown', 'unclEnUp', 'unclEnDown']:
        return "bdt_score_"+syst
    else:
        return "bdt_score_nominal"

def mllVar(syst="nominal"):
    return "dilepton_mass"

def mlljVar(syst="nominal"):
    if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown']:
        return syst+"_m_llj"
    return "nominal_m_llj"

def dRVar(syst="nominal"):
    if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown']:
        return syst+"_dR_l2j"
    else:
        return "nominal_dR_l2j"

def weight(syst="nominal"):
    return "(MET_filter > 0)*((nleadingLeptons+nsubleadingLeptons)==2)*({njets}>0)".format(njets=njetsVar(syst))

def dataWeight(sample):
    if sample=="electron":
        return "IsoElectronTrigger_flag*leadingLeptons_isElectron[0]"
    elif sample=="muon":
        return "IsoMuTrigger_flag*leadingLeptons_isMuon[0]"
    else:
        raise Exception("Unknown data name '%s'"%sample)

def mcWeight(sample, syst="nominal"):
    
    weight = str(lumiPerYear[args.year]*1000.)
    
    if syst=="puUp":
        weight += "*puweight_up"
    elif syst=="puDown":
        weight += "*puweight_down"
    else:
        weight += "*puweight_nominal"
        
    
    weight += "*(IsoElectronTrigger_flag*leadingLeptons_isElectron[0]"
    '''
    if syst=="eleEffUp":
        weight +="IsoElectronTrigger_flag*leadingLeptons_isElectron[0]*(IsoElectronTrigger_weight_trigger_up*tightElectrons_weight_reco_up*tightElectrons_weight_id_up)"
    elif syst=="eleEffDown":
        weight +="IsoElectronTrigger_flag*leadingLeptons_isElectron[0]*(IsoElectronTrigger_weight_trigger_down*tightElectrons_weight_reco_down*tightElectrons_weight_id_down)"
    else:
        weight +="IsoElectronTrigger_flag*leadingLeptons_isElectron[0]*(IsoElectronTrigger_weight_trigger_nominal*tightElectrons_weight_reco_nominal*tightElectrons_weight_id_nominal)"
    '''
    if syst=="muEffUp":
        weight +="+IsoMuTrigger_flag*leadingLeptons_isMuon[0]*(IsoMuTrigger_weight_trigger_up*tightMuons_weight_id_up*tightMuons_weight_iso_up)"
    elif syst=="muEffDown":
        weight +="+IsoMuTrigger_flag*leadingLeptons_isMuon[0]*(IsoMuTrigger_weight_trigger_down*tightMuons_weight_id_down*tightMuons_weight_iso_down)"
    else:
        weight +="+IsoMuTrigger_flag*leadingLeptons_isMuon[0]*(IsoMuTrigger_weight_trigger_nominal*tightMuons_weight_id_nominal*tightMuons_weight_iso_nominal)"
    weight+=")"
        
    return weight


def taggerScore(syst="nominal"):
    return ("({dR}>0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_Q+" \
        +"(({dR}<0.4)*(" \
            +"subleadingLeptons_isMuon[0]*hnlJet_{syst}_llpdnnx_ratio_LLP_QMU+" \
            +"subleadingLeptons_isElectron[0]*hnlJet_{syst}_llpdnnx_ratio_LLP_QE" \
        +"))").format(dR=dRVar(syst),syst=syst if syst in ['jerUp', 'jerDown', 'jesTotalUp', 'jesTotalDown'] else "nominal")


def normWeight(sample):
    xsec = -1
    integral = yields[sample]["weighted"]
    for proc in xsecs.keys():
        if sample.find(proc)>=0:
            xsec = xsecs[proc]
            break
    return "(genweight*"+str(xsec)+"/"+str(integral)+")"

def normWeightSignal(sample,coupling=1,scale=1):
    xsec = hnlXsecs[sample.replace("pt20","all")]["weights"][str(coupling)]["xsec"]["nominal"]
    if sample.find("pt20")>=0:
        xsec *= hnlFilter[sample]["weights"][str(coupling)]["eff"]

    integral = yieldsHNL[sample+"-"+str(args.year)]["LHEWeights_coupling_%i"%coupling]
    return "("+("LHEWeights_coupling_%i"%coupling)+"*genweight*"+str(xsec*scale)+"/"+str(integral)+")"

def makeHist(path,treeName,var,weight,binning,syst="nominal"):
    f = ROOT.TFile(path)
    tree = f.Get(treeName)
    #print (tree)
    if not tree:
        print ("WARNING: file '%s' does not contain a valid tree"%path)
        return None
    
    catVar =  "1*(IsoMuTrigger_flag*Leptons_muonmuon*(dilepton_charge<0))"
    catVar += "+2*(IsoElectronTrigger_flag*Leptons_electronelectron*(dilepton_charge<0))"
    catVar += "+3*(IsoMuTrigger_flag*Leptons_muonelectron*(dilepton_charge<0))"
    catVar += "+4*(IsoElectronTrigger_flag*Leptons_electronmuon*(dilepton_charge<0))"
    
    catVar += "+5*(IsoMuTrigger_flag*Leptons_muonmuon*(dilepton_charge>0))"
    catVar += "+6*(IsoElectronTrigger_flag*Leptons_electronelectron*(dilepton_charge>0))"
    catVar += "+7*(IsoMuTrigger_flag*Leptons_muonelectron*(dilepton_charge>0))"
    catVar += "+8*(IsoElectronTrigger_flag*Leptons_electronmuon*(dilepton_charge>0))"
    

    catBinning = np.linspace(-0.5,8.5,10)
    hist = ROOT.TH2F("hist_"+str(random.random()),"",len(binning)-1,binning,len(catBinning)-1,catBinning)
    hist.Sumw2()
    tree.Project(hist.GetName(),catVar+":"+var,weight)
    hist.SetDirectory(0)
    f.Close()
    return hist

def makeHistFromFolder(folder,treeName,var,weight,binning):
    files = list(filter(lambda x: x.endswith(".root"),os.listdir(folder)))
    print ("folder: ",folder,"nfiles=%i"%len(files))
    print ("var: ",var)
    print ("weight: ",weight)
    print ("binning: ",binning)

    histSum = None
    for i,f in enumerate(files):
        hist = makeHist(os.path.join(folder,f),treeName,var,weight,binning)
        if hist==None:
            continue
        if histSum==None:
            histSum = hist
        else:
            histSum.Add(hist)
        #break
    print ("Events/Integral: ",histSum.GetEntries(),"/",histSum.Integral())
    #if (histSum.Integral()/histSum.GetEntries()>1e4): #reject large MC weights
    #    histSum.Scale(0)
    print ()
    return histSum

plotCfgs ={
    "tagger_CR_boosted":{
        "var": taggerScore(args.syst),
        "cut": "({dR}<0.4)*({mll}>80.)*({njets}>0)*({nfwdjets}<1)".format(dR=dRVar(args.syst),mll=mllVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst)) ,
        "binning": np.linspace(0,1,21),
        "signals": {
            #"HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
            #    "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            #}
        },
        "blind":False
    },
    
    "tagger_CR_resolved":{
        "var": taggerScore(args.syst),
        "cut": "({dR}>0.4)*({mll}>80.)*({njets}>0)*({nfwdjets}<1)".format(dR=dRVar(args.syst),mll=mllVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst)) ,
        "binning": np.linspace(0,1,21),
        "signals": {
            #"HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
            #    "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            #}
        },
        "blind":False
    },
    
    "tagger_SR_boosted":{
        "var": taggerScore(args.syst),
        "cut": "({dR}<0.4)*({mll}>20.)*({mll}<80.)*({met}<100.)*({njets}>0)*({njets}<5)*({nfwdjets}<1)".format(dR=dRVar(args.syst),mll=mllVar(args.syst),met=metVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst)) ,
        "binning": np.linspace(0,1,21),
        "signals": {
            "HNL_dirac_pt20_ctau1p0e02_massHNL4p5_Vall1p438e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            }
        },
        "blind":False
    },
    
    "tagger_SR_resolved":{
        "var": taggerScore(args.syst),
        "cut": "({dR}>0.4)*({mll}>20.)*({mll}<80.)*({met}<100.)*({njets}>0)*({njets}<5)*({nfwdjets}<1)".format(dR=dRVar(args.syst),mll=mllVar(args.syst),met=metVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst)) ,
        "binning": np.linspace(0,1,21),
        "signals": {
            "HNL_dirac_pt20_ctau1p0e02_massHNL4p5_Vall1p438e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            }
        },
        "blind":False
    },
    
    
    "bdt_SR":{
        "var": bdtVar(args.syst),
        "cut": "({mll}>20.)*({mll}<80.)*({met}<100.)*({njets}>0)*({njets}<5)*({nfwdjets}<1)".format(mll=mllVar(args.syst),met=metVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst)) ,
        "binning": np.linspace(0,1,21),
        "signals": {
            "HNL_dirac_pt20_ctau1p0e02_massHNL4p5_Vall1p438e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            }
        },
        "blind":False
    },
    
    "mllj_SR":{
        "var": mlljVar(args.syst),
        "cut": "({mll}>20.)*({mll}<80.)*({met}<100.)*({njets}>0)*({njets}<5)*({nfwdjets}<1)".format(mll=mllVar(args.syst),met=metVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst)) ,
        "binning": np.linspace(0,300,61),
        "signals": {
            "HNL_dirac_pt20_ctau1p0e02_massHNL4p5_Vall1p438e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            }
        },
        "blind":False
    },
    
    "mll_central":{
        "var": mllVar(args.syst),
        "cut": "({mll}>10.)*({mll}<50.)*({njets}>0)*({njets}<3)*({nfwdjets}==0)*({tagger}>0.6)".format(mll=mllVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst),tagger=taggerScore(args.syst)) ,
        "binning": np.linspace(10,50,81),
        "signals": {
            "HNL_dirac_pt20_ctau1p0e02_massHNL4p5_Vall1p438e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            }
        },
        "blind":False
    },
    "mll_fwd":{
        "var": mllVar(args.syst),
        "cut": "({mll}>10.)*({mll}<50.)*({njets}>0)*({njets}<3)*({nfwdjets}==1)*({tagger}>0.6)".format(mll=mllVar(args.syst),njets=njetsVar(args.syst),nfwdjets=nfwdjetsVar(args.syst),tagger=taggerScore(args.syst)) ,
        "binning": np.linspace(10,50,81),
        "signals": {
            "HNL_dirac_pt20_ctau1p0e02_massHNL4p5_Vall1p438e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            },
            "HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03": {
                "e":2, "mu": 12, "emu": 7, "tau": 67, "all": 1
            }
        },
        "blind":False
    },
}

plotCfg = plotCfgs[args.plot]

samples = yaml.load(open("/vols/cms/mkomm/HNL/histo/config/samples.yml"))

histDict = {}

if not plotCfg['blind'] and args.syst=="nominal":
    for sampleName in ["muon","electron"]:
        histSum = None
        for run in samples[sampleName][args.year].keys():
            
            for folder in samples[sampleName][args.year][run]:
                path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21/"+str(args.year)+"/"+folder
                weight = dataWeight(sampleName)+"*"+plotCfg['cut']
                hist = makeHistFromFolder(path,"Friends",plotCfg['var'],weight,plotCfg['binning'])
                if histSum == None:
                    histSum = hist
                else:
                    histSum.Add(hist)
        histDict[sampleName] = histSum

for sampleName in ['qcd','topbkg','wjets','dyjets','vgamma']:
    histSum = None
    for proc in samples[sampleName][args.year].keys():
        for folder in samples[sampleName][args.year][proc]:
            path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21/"+str(args.year)+"/"+folder
            weight = normWeight(folder)+"*"+mcWeight(sampleName,args.syst)+"*"+plotCfg['cut']
            hist = makeHistFromFolder(path,"Friends",plotCfg['var'],weight,plotCfg['binning'])
            if histSum == None:
                histSum = hist
            else:
                histSum.Add(hist)
    histDict[sampleName] = histSum

for sampleName in plotCfg['signals']:
    folder = sampleName+"-"+str(args.year)
    path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21/"+str(args.year)+"/"+folder
    for scenarioName,coupling in plotCfg['signals'][sampleName].items():
        weight = normWeightSignal(sampleName,coupling=coupling)+"*"+mcWeight(sampleName,args.syst)+"*"+plotCfg['cut']
        hist = makeHistFromFolder(path,"Friends",plotCfg['var'],weight,plotCfg['binning'])
        histDict[sampleName+"_"+scenarioName] = hist


print ("="*60)
outputFileName = args.output+"_"+str(args.year)+"_"+args.plot+"_"+str(args.syst)+".root"
print ("write output: %s"%outputFileName)
outputFile = ROOT.TFile(outputFileName,"RECREATE")
for name,hist in histDict.items():
    hist.SetName(name)
    print ("  hist: %s"%name, "Events/Integral: ",hist.GetEntries(),"/",hist.Integral())
    hist.SetDirectory(outputFile)
outputFile.Write()
outputFile.Close()

