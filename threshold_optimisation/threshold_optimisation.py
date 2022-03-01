import ROOT
import math
import os
import pandas as pd
import ctypes
import numpy as np
import matplotlib.pyplot as plt
from cortools import mutual_information, pearson_corr
from matplotlib.colors import LogNorm
from uncertainties import ufloat
import seaborn as sns
from collections import defaultdict
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
from collections import OrderedDict as odict
import json
import pickle

bdtBinsForPlot = np.linspace(0,1,11)
taggerBinsForPlot = np.linspace(0.4,1,7)
process_name_signals = {"HNL_majorana_pt20_ctau1p0e03_massHNL2p0_Vall2p871e-03":"2GeV 1m","HNL_majorana_pt20_ctau1p0e01_massHNL8p0_Vall6p702e-04":"8GeV 10mm","HNL_majorana_pt20_ctau1p0e-01_massHNL12p0_Vall2p314e-03":"12GeV 0.1mm","HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03":"10GeV 1mm","HNL_majorana_pt20_ctau1p0e-03_massHNL16p0_Vall1p097e-02":"16GeV 0.001mm","HNL_majorana_pt20_ctau1p0e02_massHNL4p5_Vall1p016e-03":"4.5GeV 100mm"}
sigWeights = {"ee":"weight2","mumu":"weight12","emu":"weight7","etau":"weight47","mutau":"weight52"}
outputDirBaseToFormat = "abcdTaggerBdt_newInputs/thresholdOptSignif_{}"
MIN_SIGNIF_FOR_MAX = 0.0001
def discovery_score(s: float, b: float, sigmab: float) -> float:
    """
    Take s, b, sigmab and calculate the significance
    """
    if s <= 0:
        return 0.
    term_one = (s+b)*(b + sigmab ** 2)/(b ** 2 + (s+b)*sigmab ** 2)
    term_two = 1+ (s* sigmab ** 2)/(b*(b+sigmab**2))
    final = (s+b)*math.log(term_one) - b**2/sigmab**2*math.log(term_two)
    if final <= 0:
        return 0.
    else:
        return math.sqrt(2*final)




def sum_and_error(df, unityweights=0):
    if unityweights:
        sum = len(df)
        err = math.sqrt(sum)
    else:
        sum = np.sum(df.weight)
        err = math.sqrt(np.sum(df.weight ** 2))
    return max(sum, 0), err

def abcd(df, bdtCut, taggerCut, unityweights=0):
    dfAB = df[(df.bdt_score_nominal < bdtCut)]
    dfCD = df[(df.bdt_score_nominal > bdtCut)]
    dfA = dfAB[dfAB.tagger_score>taggerCut]
    dfB = dfAB[dfAB.tagger_score<taggerCut]
    dfC = dfCD[dfCD.tagger_score<taggerCut]
    dfD = dfCD[dfCD.tagger_score>taggerCut]
    sumA, errA = sum_and_error(dfA, unityweights)
    sumB, errB = sum_and_error(dfB, unityweights)
    sumC, errC = sum_and_error(dfC, unityweights)
    sumD, errD = sum_and_error(dfD, unityweights)

    if sumA == 0 or sumB == 0 or sumC == 0:
        pred = 0.
        errPred = 100.
    else:
        pred = sumA * sumC / sumB
        errPred = np.abs(pred*math.sqrt(((errA/sumA) ** 2 + (errB/sumB) ** 2 + (errC/sumC) ** 2)))
    return sumD, errD, pred, errPred

def main():

    MASS_CUT_SYM = 10
    #Min tagger score for pred
    # TAGGER_CUT_LOWER = 0.25
    optMethod = "ratio"

    cutBdtSignal = lambda df: df[(df.bdt_score_nominal>0.3)]
    cutBdtSideband = lambda df: df[(df.bdt_score_nominal<0.3)]
    cutTaggerSideband = lambda df: df[(df.tagger_score < 0.4)]
    cutTaggerSignal = lambda df: df[(df.tagger_score > 0.4)]
    cutAll = {}
    cutAll["resolved"] = lambda df: df[(df.deltamass < MASS_CUT_SYM) & (df.tagger_score>0.2)]
    cutAll["merged"] = lambda df: df[(df.deltamass < MASS_CUT_SYM) & (df.tagger_score>0.3)]
    mW = 80.3
    nMax = 3 

    # process_names_mc = [["wjets"], ["dyjets"], ["topbkg"], ["vgamma"], ["wjets", "dyjets", "topbkg", "vgamma"]]
    process_names_mc = {"mc":["wjets", "dyjets", "topbkg", "vgamma"]}#,"data":["electron","muon"]}
    process_names_data = {"data":["electron","muon"]} 
    # process_names_mc = {"data":["electron","muon"]}
    colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf',"g","b","c","r","pink"]*3
    markerStyles=['o','s', 'v', '^', '<', '>', '8', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']*3
    iC = 0
    colorDict = {}
    for process_name_signal in process_name_signals:
        for sigName,sigWeight in sigWeights.items():
            colorDict[process_name_signal,sigName] = colors[iC]
            iC += 1
    for year in ["2016","2017","2018"][:1]:

        ntuple_path = 'output_26Aug_2l/'+year
        category_names_temp = sorted(os.listdir(ntuple_path))
        signifCats = odict()
        signifCats["prompt"] = lambda df: df[(df.dxysig < 1)]
        signifCats["medium"] = lambda df: df[(df.dxysig > 1) & (df.dxysig < 10) ]
        signifCats["displaced"] = lambda df: df[(df.dxysig > 10)]
        category_signifCuts = odict()
        for category_name in category_names_temp:
            for signifCat in signifCats.keys():
                category_signifCuts[category_name+"_"+signifCat] = signifCats[signifCat]

        # BDTBinsForPlot = np.linspace(0,80,5)
        # taggerBinsForPlot = np.linspace(0.2,1,5)
        for title,process_name_set in process_names_data.items():
            totalBackgroundInBdtSideband = defaultdict(dict)
            totalBackgroundInTaggerSideband = defaultdict(dict)
            totalBackgroundInBdtSignalRegion = defaultdict(dict)
            totalBackgroundInTaggerSignalRegion = defaultdict(dict)
            totalSignalInBdtSideband = defaultdict(dict)
            totalSignalInTaggerSideband = defaultdict(dict)
            totalSignalInBdtSignalRegion = defaultdict(dict)
            totalSignalInTaggerSignalRegion = defaultdict(dict)

            coordsBestThresholds = defaultdict(dict)
            optimumOutput = defaultdict(dict)
            optimumOutputAll = defaultdict(dict)
            outputDirBase = (outputDirBaseToFormat.format(title)+"/"+year).replace(".","p")
            if not os.path.exists(outputDirBase):
                os.makedirs(outputDirBase)
            if not os.path.exists(outputDirBase+"/optima"):
                os.mkdir(outputDirBase+"/optima")
            if not os.path.exists(outputDirBase+"/optimaV2"):
                os.mkdir(outputDirBase+"/optimaV2")
            if not os.path.exists(outputDirBase+"/contamOverall"):
                os.mkdir(outputDirBase+"/contamOverall")
            for category in ["resolved", "merged"]:
                df = pd.DataFrame()
                frames = []
                for category_name in category_signifCuts.keys():
                    outputDir = outputDirBase +"/"+category+"_"+category_name
                    # if "merged_mumu_OS_displaced" not in category+"_"+category_name: continue
                    if not os.path.exists(outputDir):
                        os.makedirs(outputDir)
                    if not os.path.exists(outputDir+"/2dPlots"):
                        os.mkdir(outputDir+"/2dPlots")
                    if not os.path.exists(outputDir+"/signif"):
                        os.mkdir(outputDir+"/signif")
                    if not os.path.exists(outputDir+"/pred"):
                        os.mkdir(outputDir+"/pred")
                    if not os.path.exists(outputDir+"/contam"):
                        os.mkdir(outputDir+"/contam")

                    for process_name in process_name_set:
                        if title == "data" and category_name[0] == "e":
                            if process_name == "muon": continue
                        if title == "data" and category_name[0] == "m":
                            if process_name == "electron": continue
                        # print (category_name,title,process_name)

                        data = pd.read_pickle(os.path.join(ntuple_path, "_".join(category_name.split("_")[:-1]), process_name+".pkl"))
                        frames.append(data)

                    df = pd.concat(frames)
                    df = df[df[category]==1]
                    df= category_signifCuts[category_name](df)
                    # df = df[df.bdt_score_nominal>BDT_CUT_UNCOR]
                    df["deltamass"] = (mW-df["mass"]).abs()
                    df = cutAll[category](df)
                    plt.hist2d(df["tagger_score"],df["bdt_score_nominal"],weights=df["weight"],bins=[taggerBinsForPlot,bdtBinsForPlot])
                    plt.colorbar()
                    plt.savefig(outputDir+"/2dPlots/"+category+"_"+category_name+"_background2d.pdf")
                    plt.cla()
                    plt.clf()
                    dfAllSignals = {}
                    for process_name_signal in process_name_signals:
                        dfSig = pd.read_pickle(os.path.join(ntuple_path, "_".join(category_name.split("_")[:-1]), process_name_signal+".pkl"))
                        dfSig = dfSig[(dfSig[category]==1)]
                        dfSig["deltamass"] = (mW-dfSig["mass"]).abs()
                        dfSig = cutAll[category](dfSig)
                        dfAllSignals[process_name_signal] = dfSig
                        for sigName,sigWeight in sigWeights.items():
                            plt.hist2d(dfSig["tagger_score"],dfSig["bdt_score_nominal"],weights=dfSig[sigWeight],bins=[taggerBinsForPlot,bdtBinsForPlot])
                            plt.colorbar()
                            sigWeight = "weight12"
                            plt.savefig(outputDir+"/2dPlots/"+"_".join([category,category_name,process_name_signal,sigName])+"2d.pdf")
                            plt.cla()
                            plt.clf()
                    #Contam plots
                    dfBdtSideband = cutBdtSideband(df)
                    dfTaggerSideband = cutTaggerSideband(df)
                    dfBdtSignal = cutBdtSignal(df)
                    dfTaggerSignal = cutTaggerSignal(df)
                    dfAllSignalsBdtSideband = {process_name_signal: cutBdtSideband(dfSig) for process_name_signal,dfSig in dfAllSignals.items()}
                    dfAllSignalsTaggerSideband = {process_name_signal: cutTaggerSideband(dfSig) for process_name_signal,dfSig in dfAllSignals.items()}
                    dfAllSignalsBdtSignal = {process_name_signal: cutBdtSignal(dfSig) for process_name_signal,dfSig in dfAllSignals.items()}
                    dfAllSignalsTaggerSignal = {process_name_signal: cutTaggerSignal(dfSig) for process_name_signal,dfSig in dfAllSignals.items()}

                    taggerBackgroundDistribution,_,_ = plt.hist(dfBdtSideband["tagger_score"],bins=taggerBinsForPlot,weights=dfBdtSideband["weight"],log=True)
                    totalBackgroundInBdtSideband[category,category_name] = np.sum(dfBdtSideband["weight"])
                    taggerBackgroundDistribution[taggerBackgroundDistribution < 0] = 0
                    taggerBackgroundDistribution = np.cumsum(taggerBackgroundDistribution[::-1])[::-1]
                    for process_name_signal in process_name_signals:
                        for sigName,sigWeight in sigWeights.items():
                            # if (dfAllSignalsBdtSideband[process_name_signal][sigWeight]).sum() <= 0: continue
                            temp = plt.hist(dfAllSignalsBdtSideband[process_name_signal]["tagger_score"],bins=taggerBinsForPlot,weights=dfAllSignalsBdtSideband[process_name_signal][sigWeight],label=process_name_signals[process_name_signal]+" "+sigName,histtype="step",color = colorDict[process_name_signal,sigName])
                            totalSignalInBdtSideband[category,category_name][process_name_signal,sigName] = np.sum(dfAllSignalsBdtSideband[process_name_signal][sigWeight])
                    # plt.legend()
                    plt.savefig(outputDir+"/contam/bdtSideband_tagger"+"_".join([category,category_name])+"contam.pdf")
                    ylims = plt.gca().get_ylim()
                    plt.ylim([ylims[0],ylims[1]*20])
                    plt.cla()
                    plt.clf()

                    bdtBackgroundDistribution,_,_ = plt.hist(dfTaggerSideband["bdt_score_nominal"],bins=bdtBinsForPlot,weights=dfTaggerSideband["weight"],log=True)
                    totalBackgroundInTaggerSideband[category,category_name] = np.sum(dfTaggerSideband["weight"])
                    bdtBackgroundDistribution[bdtBackgroundDistribution < 0] = 0
                    bdtBackgroundDistribution = np.cumsum(bdtBackgroundDistribution[::-1])[::-1]
                    for process_name_signal in process_name_signals:
                        for sigName,sigWeight in sigWeights.items():
                            # if (dfAllSignalsTaggerSideband[process_name_signal][sigWeight]).sum() <= 0: continue
                            histOutput = plt.hist(dfAllSignalsTaggerSideband[process_name_signal]["bdt_score_nominal"],bins=bdtBinsForPlot,weights=dfAllSignalsTaggerSideband[process_name_signal][sigWeight],label=process_name_signals[process_name_signal]+" "+sigName,histtype="step",color=colorDict[process_name_signal,sigName])
                            totalSignalInTaggerSideband[category,category_name][process_name_signal,sigName] = np.sum(dfAllSignalsTaggerSideband[process_name_signal][sigWeight])
                    # plt.legend()
                    plt.savefig(outputDir+"/contam/taggerSideband_bdt"+"_".join([category,category_name])+"contam.pdf")
                    ylims = plt.gca().get_ylim()
                    plt.ylim([ylims[0],ylims[1]*20])
                    plt.cla()
                    plt.clf()
                    signalTaggerDistributions = {}
                    if title == "mc": 
                        plt.hist(dfBdtSignal["tagger_score"],bins=taggerBinsForPlot,weights=dfBdtSignal["weight"],log=True)
                        totalBackgroundInBdtSignalRegion[category,category_name] = np.sum(dfBdtSignal["weight"])
                    for process_name_signal in process_name_signals:
                        for sigName,sigWeight in sigWeights.items():
                            # if (dfAllSignalsBdtSignal[process_name_signal][sigWeight]).sum() <= 0: continue
                            signalTaggerDistributions[(process_name_signals[process_name_signal],sigName)],_,_ = plt.hist(dfAllSignalsBdtSignal[process_name_signal]["tagger_score"],bins=taggerBinsForPlot,weights=dfAllSignalsBdtSignal[process_name_signal][sigWeight],label=process_name_signals[process_name_signal]+" "+sigName,histtype="step")
                            signalTaggerDistributions[(process_name_signals[process_name_signal],sigName)] = np.cumsum(signalTaggerDistributions[(process_name_signals[process_name_signal],sigName)][::-1])[::-1]
                            totalSignalInBdtSignalRegion[category,category_name][process_name_signal,sigName] = np.sum(dfAllSignalsBdtSignal[process_name_signal][sigWeight])
                    # plt.legend()
                    plt.savefig(outputDir+"/contam/bdtSignal_tagger"+"_".join([category,category_name])+"contam.pdf")
                    ylims = plt.gca().get_ylim()
                    plt.ylim([ylims[0],ylims[1]*20])
                    plt.cla()
                    plt.clf()
                    signalBdtDistributions = {}

                    if title == "mc": 
                        plt.hist(dfTaggerSignal["bdt_score_nominal"],bins=bdtBinsForPlot,weights=dfTaggerSignal["weight"],log=True)
                        totalBackgroundInTaggerSignalRegion[category,category_name] = np.sum(dfTaggerSignal["weight"])
                    for process_name_signal in process_name_signals:
                        for sigName,sigWeight in sigWeights.items():
                            # if (dfAllSignalsTaggerSignal[process_name_signal][sigWeight]).sum() <= 0: continue
                            signalBdtDistributions[(process_name_signals[process_name_signal],sigName)],_,_ = plt.hist(dfAllSignalsTaggerSignal[process_name_signal]["bdt_score_nominal"],bins=bdtBinsForPlot,weights=dfAllSignalsTaggerSignal[process_name_signal][sigWeight],label=process_name_signals[process_name_signal]+" "+sigName,histtype="step")
                            signalBdtDistributions[(process_name_signals[process_name_signal],sigName)] = np.cumsum(signalBdtDistributions[(process_name_signals[process_name_signal],sigName)][::-1])[::-1]
                            totalSignalInTaggerSignalRegion[category,category_name][process_name_signal,sigName] = np.sum(dfAllSignalsTaggerSignal[process_name_signal][sigWeight])
                    # plt.legend()
                    plt.savefig(outputDir+"/contam/taggerSignal_bdt"+"_".join([category,category_name])+"contam.pdf")
                    ylims = plt.gca().get_ylim()
                    plt.ylim([ylims[0],ylims[1]*20])
                    plt.cla()
                    plt.clf()
                    #Define A/B = high BDT/low BDT, C/B = high tagger/low tagger
                    ratioCOvB = []
                    ratioAOvB = []
                    for value in bdtBackgroundDistribution:
                        total = (dfTaggerSideband["weight"]).sum()
                        if abs(value - total) < 0.01:
                            ratioAOvB.append(1)
                        else:
                            ratioAOvB.append(value/(total-value))
                    for iv,value in enumerate(taggerBackgroundDistribution):
                        total = (dfBdtSideband["weight"]).sum()
                        if abs(value - total) < 0.01:
                            ratioCOvB.append(1)
                        else:
                            ratioCOvB.append(value/(total-value))
                    iBdt = 0
                    #Do ABCD prediction
                    dfForABCD = df
                    signifScan = []
                    signifScanY = np.around(np.array([x for x in bdtBinsForPlot[:-1]]),3)
                    signifScanX = np.around(np.array([x for x in taggerBinsForPlot[:-1]]),3)
                    signifScan = defaultdict(list)
                    signalEff = defaultdict(list)
                    dPredVals = []
                    dObsVals = []
                    dPredValsUF = []
                    cPredValsUF = []
                    bPredValsUF = []
                    aPredValsUF = []
                    for bdtBinLow,bdtBinUp in zip(bdtBinsForPlot[:-1],bdtBinsForPlot[1:]):
                        iTagger = 0
                        signifScanRow = defaultdict(list)
                        signalEffRow = defaultdict(list)
                        dPredValsRow = []
                        dObsValsRow = []
                        dPredValsUFRow = []
                        cPredValsUFRow = []
                        bPredValsUFRow = []
                        aPredValsUFRow = []
                        for taggerBinLow,taggerBinUp in zip(taggerBinsForPlot[:-1],taggerBinsForPlot[1:]):
                            #ratio method
                            if optMethod == "ratio":
                                valB,errB = sum_and_error(dfForABCD[(dfForABCD.tagger_score < taggerBinLow) & (dfForABCD.bdt_score_nominal < bdtBinLow)])
                                bPred = ufloat(valB,errB)
                                if ratioCOvB[iTagger]*valB < 0:
                                    cPred = ufloat(0,100)
                                else:
                                    cPred = ufloat(ratioCOvB[iTagger]*valB,(ratioCOvB[iTagger]*valB)**0.5)
                                if ratioAOvB[iBdt]*valB < 0:
                                    aPred = ufloat(0,100)
                                else:
                                    aPred = ufloat(ratioAOvB[iBdt]*valB,(ratioAOvB[iBdt]*valB)**0.5)
                                if bPred.n <= 0:
                                    dPred = ufloat(0,100)
                                else:
                                    dPred = (aPred*cPred)/bPred
                            elif optMethod == "directABCD":
                                #direct ABCD
                                _,_,dValDirct,dErrDirect = abcd(dfForABCD,bdtBinLow,taggerBinLow)
                                dPred = ufloat(dValDirct,dErrDirect)

                            iTagger += 1
                            dPredValsUFRow.append(dPred)
                            bPredValsUFRow.append(bPred)
                            cPredValsUFRow.append(cPred)
                            aPredValsUFRow.append(aPred)
                            if title == "mc":
                                dObsValsRow.append(max((dfForABCD[(dfForABCD.tagger_score > taggerBinLow) & (dfForABCD.bdt_score_nominal > bdtBinLow)]["weight"]).sum(),0.001))
                            if dPred.n <= 0.001:
                                dPredValsRow.append(-3)
                            else:
                                dPredValsRow.append(math.log10(dPred.n))
                            for process_name_signal in process_name_signals:
                                for sigName,sigWeight in sigWeights.items():
                                    dfSig = dfAllSignals[process_name_signal]
                                    sigPred = (dfSig[(dfSig.tagger_score > taggerBinLow) & (dfSig.bdt_score_nominal > bdtBinLow)][sigWeight]).sum()
                                    sigPredA = (dfSig[(dfSig.tagger_score > taggerBinLow) & (dfSig.bdt_score_nominal < bdtBinLow)][sigWeight]).sum()
                                    sigPredC = (dfSig[(dfSig.tagger_score < taggerBinLow) & (dfSig.bdt_score_nominal > bdtBinLow)][sigWeight]).sum()
                                    aPredCorrected = ufloat(aPred.n+sigPredA,(aPred.n+sigPredA)**0.5)
                                    cPredCorrected = ufloat(cPred.n+sigPredC,(cPred.n+sigPredC)**0.5)
                                    if bPred.n < 0.001:
                                        signif = 0
                                    else:
                                        dPredCorrected = aPredCorrected*cPredCorrected/bPred
                                    if dPred.n < 0.1:# or aPredCorrected.n/aPred.n > 1.2 or cPredCorrected.n/cPred.n > 1.2:
                                        signif = 0
                                    else:
                                        s = sigPred
                                        b = dPredCorrected.n
                                        sigma = dPredCorrected.s
                                        signif = discovery_score(s,b,sigma)
                                        if signif == 0 and s > 0:
                                            print ("Warning: failure in signif calc",category,category_name,process_name_signal,sigName)
                                            print (s, b, sigma)
                                            print()
                                    signifScanRow[process_name_signal,sigName].append(signif)
                                    signalEffRow[process_name_signal,sigName].append(s/)
                        dPredVals.append(dPredValsRow)
                        dObsVals.append(dObsValsRow)
                        dPredValsUF.append(dPredValsUFRow)
                        cPredValsUF.append(cPredValsUFRow)
                        bPredValsUF.append(bPredValsUFRow)
                        aPredValsUF.append(aPredValsUFRow)
                        for process_name_signal in process_name_signals:
                            for sigName,sigWeight in sigWeights.items():
                                signifScan[process_name_signal,sigName].append(signifScanRow[process_name_signal,sigName])
                        iBdt += 1
                    dPredValsUF = np.array(dPredValsUF)
                    dPredValsUnc = []
                    for row in dPredValsUF:
                        dPredValsUncRow = []
                        for col in row:
                            dPredValsUncRow.append(col.s)
                        dPredValsUnc.append(dPredValsUncRow)
                    dPredValsUnc = np.array(dPredValsUnc)
                    cPredValsUF = np.array(cPredValsUF)
                    bPredValsUF = np.array(bPredValsUF)
                    aPredValsUF = np.array(aPredValsUF)
                    dPredVals = np.array(dPredVals)
                    dObsVals = np.array(dObsVals)
                    dPullsVals = []
                    dRatioVals = []
                    for iRow in range(len(dObsVals)):
                        dPullsValsRow = []
                        dRatioValsRow = []
                        for iCol in range(len(dObsVals[iRow])):
                            pred = dPredValsUF[iRow,iCol].n
                            unc = dPredValsUF[iRow,iCol].s
                            obs = dObsVals[iRow,iCol]
                            if pred <= 0:
                                statUnc = 1.3
                                if obs == 0:
                                    ratio = 1
                                else:
                                    ratio = -1
                            else:
                                statUnc = pred**0.5
                                ratio = obs/pred
                            if (unc**2+pred) == 0:
                                pull = 5
                            pull = (dObsVals[iRow,iCol]-pred)/(unc**2+statUnc**2)**0.5
                            dPullsValsRow.append(pull)
                            dRatioValsRow.append(ratio)
                        dPullsVals.append(dPullsValsRow)
                        dRatioVals.append(dRatioValsRow)
                            
                    ax = sns.heatmap(dPredVals,xticklabels=signifScanX,yticklabels=signifScanY,cmap="YlGnBu",vmin=-2)
                    ax.invert_yaxis()
                    plt.savefig(outputDir+"/pred"+"/"+"predForD_"+"_".join([category,category_name])+".pdf")
                    plt.cla()
                    plt.clf()
                    if title == "mc":
                        ax = sns.heatmap(np.log10(dObsVals),xticklabels=signifScanX,yticklabels=signifScanY,cmap="YlGnBu",vmin=-2)
                        ax.invert_yaxis()
                        plt.savefig(outputDir+"/pred"+"/"+"obsMCForD_"+"_".join([category,category_name])+".pdf")
                        plt.cla()
                        plt.clf()
                        ax = sns.heatmap(dRatioVals,xticklabels=signifScanX,yticklabels=signifScanY,cmap="YlGnBu",vmin=0,vmax=5)
                        ax.invert_yaxis()
                        plt.savefig(outputDir+"/pred"+"/"+"ratioPredObsForD_"+"_".join([category,category_name])+".pdf")
                        plt.cla()
                        plt.clf()
                        ax = sns.heatmap(dPullsVals,xticklabels=signifScanX,yticklabels=signifScanY,cmap="YlGnBu",vmin=-5,vmax=5)
                        ax.invert_yaxis()
                        plt.savefig(outputDir+"/pred"+"/"+"pullsObsForD_"+"_".join([category,category_name])+".pdf")
                        plt.cla()
                        plt.clf()
                    for process_name_signal in process_name_signals:
                        for sigName,sigWeight in sigWeights.items():
                            signifScan[process_name_signal,sigName] = np.array(signifScan[process_name_signal,sigName])
                            ax = sns.heatmap(signifScan[process_name_signal,sigName],xticklabels=signifScanX,yticklabels=signifScanY,cmap="YlGnBu",vmin=-3,vmax=0,annot=True,fmt=".1f")
                            ax.invert_yaxis()
                            plt.savefig(outputDir+"/signif"+"/"+"taggerSignif_"+"_".join([category,category_name,process_name_signal,sigName])+".pdf")
                            plt.cla()
                            plt.clf()
                            #Find maximal significance over threshold
                            if np.max(signifScan[process_name_signal,sigName]) > MIN_SIGNIF_FOR_MAX:
                                maximumSignifIndices = np.unravel_index(np.argsort(signifScan[process_name_signal,sigName].ravel())[-3:], signifScan[process_name_signal,sigName].shape)
                                maximumSignificances = signifScan[process_name_signal,sigName][maximumSignifIndices]
                                dataVals = dPredValsUF[maximumSignifIndices]
                                bVals = bPredValsUF[maximumSignifIndices]
                                cVals = cPredValsUF[maximumSignifIndices]
                                aVals = aPredValsUF[maximumSignifIndices]
                                coords = (signifScanX[np.array(maximumSignifIndices)[1,:]],signifScanY[np.array(maximumSignifIndices)[0,:]])
                                optimumOutput[process_name_signal,sigName][category,category_name] = (coords,maximumSignificances,dataVals,cVals,bVals,aVals)
                                #approach number 2
                                maximumSignifIndicesAll = np.unravel_index(np.argsort(signifScan[process_name_signal,sigName].ravel()), signifScan[process_name_signal,sigName].shape)
                                maximumSignificancesAll = signifScan[process_name_signal,sigName]
                                dataValsAll = dPredValsUF
                                bValsAll = bPredValsUF
                                cValsAll = cPredValsUF
                                aValsAll = aPredValsUF
                                coordsAll = (signifScanX,signifScanY)
                                optimumOutputAll[process_name_signal,sigName][category,category_name] = (coordsAll,maximumSignificancesAll,dataValsAll,cValsAll,bValsAll,aValsAll)
                            else:
                                optimumOutputAll[process_name_signal,sigName][category,category_name] = ((signifScanX,signifScanY),signifScan[process_name_signal,sigName],dPredValsUF,cPredValsUF,bPredValsUF,aPredValsUF)

        yLabels = []
        xLabels = []
        contam2DBdt = []
        contam2DTagger = []
        if title == "mc":
            for category,category_name in totalSignalInTaggerSignalRegion:
                yLabels.append(category+"\n"+category_name)
                contam2DBdtRow = []
                contam2DTaggerRow = []
                for process_name_signal,sigName in totalSignalInTaggerSignalRegion[category,category_name]:
                    if totalSignalInBdtSignalRegion[category,category_name][process_name_signal,sigName] < 0.1:
                        doubleRatioBdt = 1E-5
                    else:
                        doubleRatioBdt = (totalSignalInBdtSideband[category,category_name][process_name_signal,sigName]/totalBackgroundInTaggerSideband[category,category_name])/(totalSignalInBdtSignalRegion[category,category_name][process_name_signal,sigName]/totalBackgroundInTaggerSignalRegion[category,category_name])
                    if totalSignalInTaggerSignalRegion[category,category_name][process_name_signal,sigName] < 0.1:
                        doubleRatioTagger = 1E-5
                    else:
                        doubleRatioTagger = (totalSignalInTaggerSideband[category,category_name][process_name_signal,sigName]/totalBackgroundInTaggerSideband[category,category_name])/(totalSignalInTaggerSignalRegion[category,category_name][process_name_signal,sigName]/totalBackgroundInTaggerSignalRegion[category,category_name])
                    doubleRatioBdt = max(doubleRatioBdt,1E-5)
                    doubleRatioTagger = max(doubleRatioTagger,1E-5)
                    contam2DBdtRow.append(doubleRatioBdt)
                    contam2DTaggerRow.append(doubleRatioTagger)
                contam2DBdt.append(contam2DBdtRow)
                contam2DTagger.append(contam2DTaggerRow)
            contam2DBdt = np.array(contam2DBdt)
            contam2DTagger = np.array(contam2DTagger)
            for process_name_signal,sigName in totalSignalInTaggerSignalRegion[category,category_name]:
                xLabels.append(process_name_signals[process_name_signal]+"\n"+sigName)
            sns.set(font_scale=0.4)
            ax = sns.heatmap(np.log10(contam2DBdt),xticklabels=xLabels,yticklabels=yLabels,cmap="YlGnBu",vmin=-5,vmax=1)
            ax.invert_yaxis()
            plt.savefig(outputDirBase+"/contamOverall"+"/"+"bdtSideband.pdf")
            plt.cla()
            plt.clf()
            ax = sns.heatmap(np.log10(contam2DTagger),xticklabels=xLabels,yticklabels=yLabels,cmap="YlGnBu",vmin=-5,vmax=1)
            ax.invert_yaxis()
            plt.savefig(outputDirBase+"/contamOverall"+"/"+"taggerSideband.pdf")
            plt.cla()
            plt.clf()
            sns.set(font_scale=1.0)
        skipCat = defaultdict(dict)
        totalSignifs = {}
        plt.rcParams["figure.figsize"] = (5,8)
        pickle.dump(optimumOutputAll,open(outputDirBase+"/optimumOutputAll.pkl","wb"))
        for process_name_signal,sigName in optimumOutput:
            #which categories to consider
            totalSignif = 0
            norm = 0
            for cat in optimumOutput[process_name_signal,sigName]:
                totalSignif += optimumOutput[process_name_signal,sigName][cat][1][-1]
                norm += 1
            totalSignif = totalSignif#/norm**0.5
            totalSignifs[process_name_signal,sigName] = totalSignif
            for cat in optimumOutput[process_name_signal,sigName]:
                if optimumOutput[process_name_signal,sigName][cat][1][-1] < totalSignif*0.01:
                    skipCat[cat][process_name_signal,sigName] = True
                else:
                    skipCat[cat][process_name_signal,sigName] = False
        for cat in skipCat:
            if all(value for value in skipCat[cat].values()):
                for proc in skipCat[cat]:
                    skipCat[cat][proc] = False
        for cat in skipCat:
            maxForHistTotal = np.zeros((len(signifScanY),len(signifScanX)))
            maxForHistTotalWeighted = np.zeros((len(signifScanY),len(signifScanX)))
            norm = np.sum([value==False for value in skipCat[cat].values()])
            for (process_name_signal,sigName) in skipCat[cat]:
                if skipCat[cat][process_name_signal,sigName]:
                    continue

                coords = optimumOutput[process_name_signal,sigName][cat][0]
                maximumSignificances = optimumOutput[process_name_signal,sigName][cat][1]
                dataVals = optimumOutput[process_name_signal,sigName][cat][2]
                cVals = optimumOutput[process_name_signal,sigName][cat][3]
                bVals = optimumOutput[process_name_signal,sigName][cat][4]
                aVals = optimumOutput[process_name_signal,sigName][cat][5]
                for i in range(nMax):
                    if maximumSignificances[i] < MIN_SIGNIF_FOR_MAX: continue
                    if title == "mc":
                        plt.scatter(coords[0][i],coords[1][i],marker=markerStyles[i],color=colorDict[process_name_signal,sigName],label="{0} {1} max {2} ({3} of total) at {4},{5} (bkg: {6}*{7}/{8} = ${9}\pm{10}$)".format(process_name_signals[process_name_signal],sigName,"$%.2f\sigma$"%maximumSignificances[i],"%.2f"%(maximumSignificances[i]/totalSignifs[process_name_signal,sigName]),coords[0][i],coords[1][i],"%.1f"%aVals[i].n,"%.1f"%cVals[i].n,"%.1f"%bVals[i].n,"%.2f"%(dataVals[i].n),"%.2f"%(dataVals[i].s)),alpha=0.4)
                    else:
                        plt.scatter(coords[0][i],coords[1][i],marker=markerStyles[i],color=colorDict[process_name_signal,sigName],label="{0} {1} max {2} ({3} of total) at {4},{5} (bkg: ${9}\pm{10}$)".format(process_name_signals[process_name_signal],sigName,"$%.2f\sigma$"%maximumSignificances[i],"%.2f"%(maximumSignificances[i]/totalSignifs[process_name_signal,sigName]),coords[0][i],coords[1][i],"%.1f"%aVals[i].n,"%.1f"%cVals[i].n,"??","%.2f"%(dataVals[i].n),"%.2f"%(dataVals[i].s)),alpha=0.4)
                maximumSignificancesAll = optimumOutputAll[process_name_signal,sigName][cat][1]
                for ix,iy in np.ndindex(maxForHistTotal.shape):
                    maxSig = maximumSignificancesAll[ix,iy]
                    if maxSig > np.max(maximumSignificancesAll)*0.7:
                        maxForHistTotal[ix,iy] += 1./norm
                    maxForHistTotalWeighted[ix,iy] += (maxSig/totalSignifs[process_name_signal,sigName])/norm

                
            plt.xlim([signifScanX[0],signifScanX[-1]])
            plt.ylim([signifScanY[0],signifScanY[-1]])
            plt.gca().invert_yaxis()
            plt.legend(prop={'size': 5},loc="upper center")
            plt.savefig(outputDirBase+"/optima"+"/"+"optima_"+"_".join(cat)+".pdf")
            plt.cla()
            plt.clf()
            ax = sns.heatmap(maxForHistTotal,xticklabels=signifScanX,yticklabels=signifScanY,vmin=0,cmap="YlGnBu")
            ax.invert_yaxis()
            plt.savefig(outputDirBase+"/optimaV2"+"/"+"optimaV2_"+"_".join(cat)+".pdf")
            plt.cla()
            plt.clf()
            maximumWeightIndices = np.unravel_index(np.argsort(maxForHistTotalWeighted.ravel())[-1:], maxForHistTotalWeighted.shape)
            dataVals = optimumOutputAll[process_name_signal,sigName][cat][2]
            dVal = dataVals[maximumWeightIndices[0],[maximumWeightIndices[1]]][0][0]
            cVals = optimumOutputAll[process_name_signal,sigName][cat][3]
            cVal = cVals[maximumWeightIndices[0],[maximumWeightIndices[1]]][0][0]
            bVals = optimumOutputAll[process_name_signal,sigName][cat][4]
            bVal = bVals[maximumWeightIndices[0],[maximumWeightIndices[1]]][0][0]
            aVals = optimumOutputAll[process_name_signal,sigName][cat][5]
            aVal = aVals[maximumWeightIndices[0],[maximumWeightIndices[1]]][0][0]
            abcdBestThresholds = ("\nbkg at {5},{6}: {0}*{1}/{2} = ${3}\pm{4}$").format("%.1f"%aVal.n,"%.1f"%cVal.n,"%.1f"%bVal.n,"%.2f"%(dVal.n),"%.2f"%(dVal.s),"%.2f"%signifScanX[np.array(maximumWeightIndices)[1,:]][0],"%.1f"%signifScanY[np.array(maximumWeightIndices)[0,:]][0])
            ax = sns.heatmap(maxForHistTotalWeighted,xticklabels=signifScanX,yticklabels=signifScanY,vmin=0,cmap="YlGnBu",annot=True,fmt=".1f")
            ax.invert_yaxis()
            plt.title(" ".join(cat).replace("_"," ")+abcdBestThresholds)
            plt.text(0.4,17.5,abcdBestThresholds)
            plt.savefig(outputDirBase+"/optimaV2"+"/"+"optimaV2Weighted_"+"_".join(cat)+".pdf")
            plt.cla()
            plt.clf()
            coordsBestThresholds[cat[1]][cat[0]] = (signifScanX[np.array(maximumWeightIndices)[1,:]][0],signifScanY[np.array(maximumWeightIndices)[0,:]][0])
        with open(outputDirBase+"/coordsBestThresholds_{0}.json".format(year),"w") as f: 
            json.dump(coordsBestThresholds,f)


if __name__ == "__main__":
    main()
