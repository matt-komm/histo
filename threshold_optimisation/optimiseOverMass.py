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
from collections import OrderedDict as odict
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import json
import pickle
from threshold_optimisation import bdtBinsForPlot,taggerBinsForPlot,process_name_signals,sigWeights,MIN_SIGNIF_FOR_MAX

# inputDirectories = {0.4:"newBenchmarks/thresholdOptVSignif_mc_0p4/",0.55:"newBenchmarks/thresholdOptVSignif_mc_0p55",0.7:"newBenchmarks/thresholdOptVSignif_mc_0p7/"}
inputDirectories = {10.:"abcdTaggerBdt_newInputs{}/thresholdOptSignif_{}/",11.:"abcdTaggerBdt_newInputs{}/thresholdOptSignif_{}/"}
# inputDirectories = {0.4:"inputBugFixLowerMassOnly{}/thresholdOptSignif_{}_0p4/",0.55:"inputBugFixLowerMassOnly{}/thresholdOptSignif_{}_0p55",0.7:"inputBugFixLowerMassOnly{}/thresholdOptSignif_{}_0p7"}
# signifScanZ = [0.4,0.8]
tag = "AllData"
dType = "data"
outputDir = "summary_newVars/{}/outputSummaryVSignif".format(tag)
if not os.path.exists(outputDir+"/overBdt"):
    os.makedirs(outputDir+"/overBdt")
finalOutputTaggerDict = defaultdict(dict)
finalOutputBdtDict = defaultdict(dict)
for year in ["2016","2017","2018"]:
    bestThresholdsPerBdt = {}
    optimumOutputAll = {}
    ntuple_path = 'output_28May21_2l/'+year
    category_names_temp = sorted(os.listdir(ntuple_path))
    category_names = []
    for category_name in category_names_temp:
        for signifCat in ["prompt","medium","displaced"]:
            category_names.append(category_name+"_"+signifCat)
    bestThresholdsOverall = defaultdict(dict)
    for bdtCut,inputDirectory in inputDirectories.items():
        xLabels = []
        yLabels = []
        bestThresholdsPerBdt[bdtCut] = json.load(open(inputDirectory.format(tag,dType)+"/"+year+"/"+"coordsBestThresholds_{0}.json".format(year),'r'))
        optimumOutputAll[bdtCut] = pickle.load(open(inputDirectory.format(tag,dType)+"/"+year+"/"+"optimumOutputAll.pkl",'rb'))
        overallSignifcancesPerCategory = []
        for category_name in category_names:
            for category in ["merged","resolved"]:
                yLabels.append("\n".join([category,category_name]))
        maxSignificanceSummaryPlot = []
        bdtRealSummaryPlot = []
        taggerCutsSummaryPlot = []
        for process_name_signal in process_name_signals:
            for sigName in sigWeights:
                proc = (process_name_signal,sigName)
                xLabels.append("\n".join([process_name_signals[proc[0]],proc[1]]))
                maxSignificanceSummaryPlotRow = []
                taggerCutsSummaryPlotRow = []
                bdtRealSummaryPlotRow = []
                for category_name in category_names:
                    for category in ["merged","resolved"]:
                        if proc not in optimumOutputAll[bdtCut].keys():
                            maxSignificanceSummaryPlotRow.append(0)
                        else:
                            if category_name in bestThresholdsPerBdt[bdtCut]:
                                print (category_name)
                                if category not in bestThresholdsPerBdt[bdtCut][category_name]:
                                    maxSignificanceSummaryPlotRow.append(0)
                                    taggerCutsSummaryPlotRow.append(-1)
                                    bdtRealSummaryPlotRow.append(-1)
                                else:
                                    vals = bestThresholdsPerBdt[bdtCut][category_name][category]
                                    signifScanX,signifScanY = optimumOutputAll[bdtCut][proc][category,category_name][0] 
                                    indices = [list(signifScanX).index(vals[0]),list(signifScanY).index(vals[1])]
                                    significances = optimumOutputAll[bdtCut][proc][category,category_name][1]
                                    maxSignificanceSummaryPlotRow.append(significances[indices[1],indices[0]])
                                    taggerCutsSummaryPlotRow.append(vals[0])
                                    bdtRealSummaryPlotRow.append(vals[1])
                            else:
                                maxSignificanceSummaryPlotRow.append(0)
                                taggerCutsSummaryPlotRow.append(-1)
                                bdtRealSummaryPlotRow.append(-1)
                maxSignificanceSummaryPlot.append(maxSignificanceSummaryPlotRow)
                taggerCutsSummaryPlot.append(taggerCutsSummaryPlotRow)
                bdtRealSummaryPlot.append(bdtRealSummaryPlotRow)
        bdtRealSummaryPlot = np.array(bdtRealSummaryPlot)
        taggerCutsSummaryPlot = np.array(taggerCutsSummaryPlot)
        maxSignificanceSummaryPlot = np.array(maxSignificanceSummaryPlot)
        maxSignificanceSummaryPlot[maxSignificanceSummaryPlot < 0] = 0 
        sns.set(font_scale=0.3)
        sns.heatmap(maxSignificanceSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".2f")
        plt.savefig(outputDir+"/summarySignificance_BDT_{0}_{1}.pdf".format(bdtCut,year))
        plt.cla()
        plt.clf()
        sns.heatmap(bdtRealSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
        plt.savefig(outputDir+"/bdtReal_BDT_{0}_{1}.pdf".format(bdtCut,year))
        plt.cla()
        plt.clf()
        sns.heatmap(taggerCutsSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
        plt.savefig(outputDir+"/taggerCuts_BDT_{0}_{1}.pdf".format(bdtCut,year))
        plt.cla()
        plt.clf()

    optimumOutputOverBdt = defaultdict(dict)
    for process_name_signal in process_name_signals.keys():
        for sigName in sigWeights:
            for category_name in category_names:
                for category in ["merged","resolved"]:
                    proc = (process_name_signal,sigName)
                    optimumOutputOverBdtTemp = []
                    cat = (category,category_name)
                    for i in range(6): optimumOutputOverBdtTemp.append([])
                    for bdtCut,inputDirectory in inputDirectories.items():
                        if len(optimumOutputOverBdtTemp[0]) == 0:
                            optimumOutputOverBdtTemp[0] = list(optimumOutputAll[bdtCut][proc][cat][0])
                            optimumOutputOverBdtTemp[0].append([bdtCut])
                            for i in range(1,len(optimumOutputAll[bdtCut][proc][cat])):
                                optimumOutputOverBdtTemp[i] = optimumOutputAll[bdtCut][proc][cat][i]
                        else:
                            optimumOutputOverBdtTemp[0] = list(optimumOutputOverBdtTemp[0])
                            optimumOutputOverBdtTemp[0][-1].append(bdtCut)
                            if len(optimumOutputOverBdtTemp[0][-1]) == 2:
                                for i in range(1,len(optimumOutputAll[bdtCut][proc][cat])):
                                    optimumOutputOverBdtTemp[i] = np.stack([optimumOutputOverBdtTemp[i],optimumOutputAll[bdtCut][proc][cat][i]],-1)
                            else:
                                for i in range(1,len(optimumOutputAll[bdtCut][proc][cat])):
                                    optimumOutputOverBdtTemp[i] = np.dstack((optimumOutputOverBdtTemp[i],optimumOutputAll[bdtCut][proc][cat][i]))
                                    print (optimumOutputOverBdtTemp[i].shape)

                    optimumOutputOverBdtTemp[0][-1] = np.array(optimumOutputOverBdtTemp[0][-1])
                    optimumOutputOverBdtTemp[0] = tuple(optimumOutputOverBdtTemp[0])
                    optimumOutputOverBdt[proc][cat] = optimumOutputOverBdtTemp
    skipCat = defaultdict(dict)
    totalSignif = defaultdict(float)
    norm = defaultdict(float)
    signifScanX,signifScanY,signifScanZ = optimumOutputOverBdt[proc][cat][0] 

    for proc in optimumOutputOverBdt:
        for cat in optimumOutputOverBdt[proc]:
            if np.max(optimumOutputOverBdt[proc][cat][1]) < MIN_SIGNIF_FOR_MAX:
                norm[cat] += 0
                totalSignif[proc] += 0
            else:
                norm[cat] += 1
                totalSignif[proc] += np.max(optimumOutputOverBdt[proc][cat][1])
    #Norm is over all signal processes that contribue to the bin
    for cat in norm:
        maxForHistTotalWeighted = np.zeros((len(signifScanY),len(signifScanX),len(signifScanZ)))
        if norm[cat] > 0:
            for proc in totalSignif.keys():
                if totalSignif[proc] == 0: continue
                if np.max(optimumOutputOverBdt[proc][cat][1]) < MIN_SIGNIF_FOR_MAX: continue
                maximumSignificancesAll = optimumOutputOverBdt[proc][cat][1]
                print (maximumSignificancesAll.shape)
                for ix,iy,iz in np.ndindex(maxForHistTotalWeighted.shape):
                    maxSig = maximumSignificancesAll[ix,iy,iz]
                    maxForHistTotalWeighted[ix,iy,iz] += (maxSig/totalSignif[proc])/norm[cat]
        maximumWeightIndices = np.unravel_index(np.argsort(maxForHistTotalWeighted.ravel())[-1:], maxForHistTotalWeighted.shape)
        dataVals = optimumOutputOverBdt[proc][cat][2]
        dVal = dataVals[maximumWeightIndices[0],maximumWeightIndices[1],maximumWeightIndices[2]][0]
        cVals = optimumOutputOverBdt[proc][cat][3]
        cVal = cVals[maximumWeightIndices[0],maximumWeightIndices[1],maximumWeightIndices[2]][0]
        bVals = optimumOutputOverBdt[proc][cat][4]
        bVal = bVals[maximumWeightIndices[0],maximumWeightIndices[1],maximumWeightIndices[2]][0]
        aVals = optimumOutputOverBdt[proc][cat][5]
        aVal = aVals[maximumWeightIndices[0],maximumWeightIndices[1],maximumWeightIndices[2]][0]
        abcdBestThresholds = ("\nbkg at {5},{6}: {0}*{1}/{2} = ${3}\pm{4}$").format("%.1f"%aVal.n,"%.1f"%cVal.n,"%.1f"%bVal.n,"%.2f"%(dVal.n),"%.2f"%(dVal.s),"%.2f"%signifScanX[np.array(maximumWeightIndices)[1,:]][0],"%.1f"%signifScanY[np.array(maximumWeightIndices)[0,:]][0])
        bestThresholdsOverall[cat[1]][cat[0]] = (signifScanX[np.array(maximumWeightIndices)[1,:]][0],signifScanY[np.array(maximumWeightIndices)[0,:]][0],signifScanZ[np.array(maximumWeightIndices)[2,:]][0])
    maxSignificanceOverBdtSummaryPlot = []
    taggerCutsOverBdtSummaryPlot = []
    bdtRealOverBdtSummaryPlot = []
    bdtCutsOverBdtSummaryPlot = []
    for process_name_signal in process_name_signals:
        for sigName in sigWeights:
            proc = (process_name_signal,sigName)
            maxSignificanceOverBdtSummaryPlotRow = []
            taggerCutsOverBdtSummaryPlotRow = []
            bdtRealOverBdtSummaryPlotRow = []
            bdtCutsOverBdtSummaryPlotRow = []
            for category_name in category_names:
                for category in ["merged","resolved"]:
                    vals = bestThresholdsOverall[category_name][category]
                    indices = [list(signifScanX).index(vals[0]),list(signifScanY).index(vals[1]),list(signifScanZ).index(vals[2])]
                    significances = optimumOutputOverBdt[proc][category,category_name][1]
                    maxSignificanceOverBdtSummaryPlotRow.append(significances[indices[1],indices[0],indices[2]])
                    taggerCutsOverBdtSummaryPlotRow.append(vals[0])
                    bdtRealOverBdtSummaryPlotRow.append(vals[1])
                    bdtCutsOverBdtSummaryPlotRow.append(vals[2])
                    disp = category_name.split("_")[-1]
                    finalOutputBdtDict[year+" "+category+"\n"+" ".join(category_name.split("_")[:-1])][disp] = bdtRealOverBdtSummaryPlotRow[-1]
                    finalOutputTaggerDict[year+" "+category+"\n"+" ".join(category_name.split("_")[:-1])][disp] = taggerCutsOverBdtSummaryPlotRow[-1]
            maxSignificanceOverBdtSummaryPlot.append(maxSignificanceOverBdtSummaryPlotRow)
            taggerCutsOverBdtSummaryPlot.append(taggerCutsOverBdtSummaryPlotRow)
            bdtRealOverBdtSummaryPlot.append(bdtRealOverBdtSummaryPlotRow)
            bdtCutsOverBdtSummaryPlot.append(bdtCutsOverBdtSummaryPlotRow)
    bdtRealOverBdtSummaryPlot = np.array(bdtRealOverBdtSummaryPlot)
    taggerCutsOverBdtSummaryPlot = np.array(taggerCutsOverBdtSummaryPlot)
    bdtCutsOverBdtSummaryPlot = np.array(bdtCutsOverBdtSummaryPlot)
    maxSignificanceOverBdtSummaryPlot = np.array(maxSignificanceOverBdtSummaryPlot)
    maxSignificanceOverBdtSummaryPlot[maxSignificanceOverBdtSummaryPlot < 0] = 0 
    sns.set(font_scale=0.3)
    sns.heatmap(maxSignificanceOverBdtSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
    plt.savefig(outputDir+"/overBdt/summarySignificanceOverBDT_{}.pdf".format(year))
    plt.cla()
    plt.clf()
    sns.heatmap(bdtRealOverBdtSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
    plt.savefig(outputDir+"/overBdt/bdtRealOverBDT_{}.pdf".format(year))
    plt.cla()
    plt.clf()
    sns.heatmap(taggerCutsOverBdtSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
    plt.savefig(outputDir+"/overBdt/taggerCutsOverBDT_{}.pdf".format(year))
    plt.cla()
    plt.clf()
    sns.heatmap(bdtCutsOverBdtSummaryPlot,xticklabels=yLabels,yticklabels=xLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
    plt.savefig(outputDir+"/overBdt/bdtCutsOverBDT_{}.pdf".format(year))
    plt.cla()
    plt.clf()
    with open(outputDir+"/coordsBestThresholdsIncBdt_{0}.json".format(year),"w") as f: 
        json.dump(bestThresholdsOverall,f)
finalOutputBdt = []
finalOutputTagger = []
yLabels = odict()
xLabels = []
for year in ["2016","2017","2018"]:
    for category in ["resolved","merged"]:
        for lep in ["ee","mumu","mue","emu"]:
            for sign in ["OS","SS"]:
                yLabels[year+" "+category+"\n"+lep+ " "+sign] = year +" "+ category+"\n"+lep+" "+sign
for xLabel in ["prompt","medium","displaced"]:
    xLabels.append(xLabel)

for keyY,yLabel in yLabels.items():
    finalOutputBdtRow = []
    finalOutputTaggerRow = []
    for xLabel in xLabels:
        finalOutputBdtRow.append(finalOutputBdtDict[yLabel][xLabel])
        finalOutputTaggerRow.append(finalOutputTaggerDict[yLabel][xLabel])
    finalOutputBdt.append(finalOutputBdtRow)
    finalOutputTagger.append(finalOutputTaggerRow)

sns.heatmap(finalOutputBdt,xticklabels=xLabels,yticklabels=list(yLabels.keys()),cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
plt.savefig(outputDir+"/finalOutputBDT.pdf")
plt.cla()
plt.clf()


sns.heatmap(finalOutputTagger,xticklabels=xLabels,yticklabels=yLabels,cmap="YlGnBu",vmin=0,vmax=1,annot=True,fmt=".1f")
plt.savefig(outputDir+"/finalOutputTagger.pdf")
plt.cla()
plt.clf()
