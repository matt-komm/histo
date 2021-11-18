import ROOT
import math
import random
import numpy as np
from histo import style

ROOT.gROOT.SetBatch(True)

colors = []
    
def newColorRGB(red,green,blue):
    color=ROOT.TColor(red,green,blue,1.0)
    colors.append(color)
    return color

newColorRGB.colorindex = 301

mcStyles = {
    'topbkg': {'legend': "t#bar{t} / t", 'fill': newColorRGB(1.,0.8,0.0)},
    'wjets': {'legend': "W+jets", 'fill': newColorRGB(0.33,0.75,0.35)},
    'dyjets':{'legend': "Z/#gamma*+jets", 'fill': newColorRGB(0.3,0.75,0.95)},
    'vgamma': {'legend': "V#gamma*+jets", 'fill': newColorRGB(0.73,0.25,0.96)},
    'qcd': {'legend': "Multijet", 'fill': newColorRGB(0.85,0.85,0.85)}
}

for name,mcStyle in mcStyles.items():
    fillColor = mcStyle['fill']
    mcStyle['line'] = newColorRGB(
        fillColor.GetRed()*0.6,
        fillColor.GetGreen()*0.6,
        fillColor.GetBlue()*0.6,
    )


class Plot():
    def __init__(self,
        plot,
        title,
        combine = lambda h: h.ProjectionX(h.GetName()+"tot"),
        extraTitle="",
        logy=False,
        unit="",
        binRange = [0,1],
        outputSuffix = "",
        procs = ['topbkg','wjets','dyjets','vgamma','qcd'],
        path='/vols/cms/mkomm/HNL/histo/plot_paper/hists/hist'
    ):
        self.plot = plot
        self.title = title
        self.combine = combine
        self.logy = logy
        self.unit = unit
        self.binRange = binRange
        self.outputSuffix = outputSuffix
        self.procs = procs
        self.path = path
        self.signals = []

    def addSignal(self, name, scale, legend, style=[2,2,ROOT.kRed+1]):
        self.signals.append({
            'name':name,
            'scale': scale,
            'legend': legend,
            'style': style
        })
        
    def getMCSignal(self,signal,syst='nominal'):
        mcHistSum = None
        for year in [2016,2017,2018]:
            filepath = self.path+"_"+str(year)+"_"+self.plot+"_"+syst+".root"
            f = ROOT.TFile(filepath)
            if f==None:
                raise Exception("Cannot open file '%s'"%filepath)
                
            hist = f.Get(signal['name'])
            hist = self.combine(hist)
            hist.SetDirectory(0)
            hist.SetFillStyle(0)
            hist.SetLineWidth(signal['style'][0])
            hist.SetLineStyle(signal['style'][1])
            hist.SetLineColor(signal['style'][2])
    
            if mcHistSum==None:
                mcHistSum = hist.Clone(self.plot+str(random.random())+"sum")
                mcHistSum.SetDirectory(0)
            else:
                mcHistSum.Add(hist)
            f.Close()
        mcHistSum.Scale(signal['scale'])
        return mcHistSum
        
    def getMCStackSumDict(self,syst='nominal'):
        mcHistDict = {}
        mcHistSum = None
        for year in [2016,2017,2018]:
            filepath = self.path+"_"+str(year)+"_"+self.plot+"_"+syst+".root"
            f = ROOT.TFile(filepath)
            if f==None:
                raise Exception("Cannot open file '%s'"%filepath)
                
            for i,sample in enumerate(self.procs):
                hist = f.Get(sample)
                hist = self.combine(hist)
                hist.SetDirectory(0)
                hist.SetFillColor(mcStyles[sample]['fill'].GetNumber())
                hist.SetLineColor(mcStyles[sample]['line'].GetNumber())
                hist.SetLineWidth(2)
                if sample not in mcHistDict.keys():
                    mcHistDict[sample] = hist.Clone(self.plot+sample+str(random.random()))
                    mcHistDict[sample].SetDirectory(0)
                else:
                    mcHistDict[sample].Add(hist)
                if mcHistSum==None:
                    mcHistSum = hist.Clone(self.plot+str(random.random())+"sum")
                    mcHistSum.SetDirectory(0)
                else:
                    mcHistSum.Add(hist)
            f.Close()
        
        mcStack = ROOT.THStack()
        for sample in self.procs:
            mcStack.Add(mcHistDict[sample])
            
        return mcStack,mcHistDict,mcHistSum
        
    def getMC(self):
        mcStackNominal, mcHistDictNominal, mcHistSumNominal = self.getMCStackSumDict()
        
        totalUnc2 = np.zeros(mcHistSumNominal.GetNbinsX())
        for syst in ['jesTotal','jer','unclEn','pu','muEff']:
            _, _, mcHistSumSystUp = self.getMCStackSumDict(syst+"Up")
            _, _, mcHistSumSystDown = self.getMCStackSumDict(syst+"Down")
            for ibin in range(mcHistSumNominal.GetNbinsX()):
                nom = mcHistSumNominal.GetBinContent(ibin+1)
                upDiff = math.fabs(mcHistSumSystUp.GetBinContent(ibin+1)-nom)
                downDiff = math.fabs(mcHistSumSystDown.GetBinContent(ibin+1)-nom)
                diff = 0.5*(upDiff+downDiff) #take average
                totalUnc2[ibin]+=diff**2
                
        for ibin in range(mcHistSumNominal.GetNbinsX()):
            mcHistSumNominal.SetBinError(ibin+1,math.sqrt(mcHistSumNominal.GetBinError(ibin+1)**2+totalUnc2[ibin])) #add mcstat err
                
        return mcStackNominal, mcHistDictNominal, mcHistSumNominal

    def __call__(self):
        
        cv = style.makeCanvas("cv"+str(random.random()),700,670)
        cv.Divide(1,2,0,0)
        cv.GetPad(1).SetPad(0.0, 0.0, 1.0, 1.0)
        cv.GetPad(1).SetFillStyle(4000)
        cv.GetPad(2).SetPad(0.0, 0.00, 1.0,1.0)
        cv.GetPad(2).SetFillStyle(4000)

        cvxmin=0.13
        cvxmax=0.97
        cvymin=0.145
        cvymax=0.92
        resHeight=0.35

        for i in range(1,3):
            #for the canvas:
            cv.GetPad(i).SetBorderMode(0)
            cv.GetPad(i).SetGridx(False)
            cv.GetPad(i).SetGridy(False)


            #For the frame:
            cv.GetPad(i).SetFrameBorderMode(0)
            cv.GetPad(i).SetFrameBorderSize(1)
            cv.GetPad(i).SetFrameFillColor(0)
            cv.GetPad(i).SetFrameFillStyle(0)
            cv.GetPad(i).SetFrameLineColor(1)
            cv.GetPad(i).SetFrameLineStyle(1)
            cv.GetPad(i).SetFrameLineWidth(1)

            # Margins:
            cv.GetPad(i).SetLeftMargin(cvxmin)
            cv.GetPad(i).SetRightMargin(1-cvxmax)

            # For the Global title:
            cv.GetPad(i).SetTitle("")

            # For the axis:
            cv.GetPad(i).SetTickx(1)  # To get tick marks on the opposite side of the frame
            cv.GetPad(i).SetTicky(1)

            # Change for log plots:
            cv.GetPad(i).SetLogx(0)
            cv.GetPad(i).SetLogy(0)
            cv.GetPad(i).SetLogz(0)

        cv.GetPad(2).SetTopMargin(1-cvymax)
        cv.GetPad(2).SetBottomMargin(resHeight)
        cv.GetPad(1).SetTopMargin(1-resHeight)
        cv.GetPad(1).SetBottomMargin(cvymin)

        cv.cd(2)
        
        
        mcStackNominal, mcHistDictNominal, mcHistSumNominal = self.getMC()
        
        if self.unit!="":
            xaxisTitle = self.title+" ("+self.unit+")"
            
            binwidth = mcHistSumNominal.GetBinWidth(1)
            binexp = math.floor(math.log10(binwidth))
            if binexp>3 or binexp<-3:
                yaxisTitle = "Events / %4.3f#upoint 10#scale[0.7]{#lower[-0.7]{%i}} %s"%(binwidth/(10**binexp),binexp,unit)
            elif binexp>1:
                yaxisTitle = "Events / %3f %s"%(binwidth,unit)
            else:
                yaxisTitle = "Events / %4.2f %s"%(binwidth,unit)
        else:
            xaxisTitle = self.title
            yaxisTitle = "Events / bins"
        
        ymax = mcHistSumNominal.GetMaximum()

        if self.logy:
            axis=ROOT.TH2F(
                "axis"+str(random.randint(0,99999)),";;"+yaxisTitle,
                50,self.binRange[0],self.binRange[1],
                50,0.07,10**((1.63+0.2*len(self.signals))*math.log10(max([ymax,1.0])))
            )
        else:
            axis=ROOT.TH2F(
                "axis"+str(random.randint(0,99999)),";;"+yaxisTitle,
                50,self.binRange[0],self.binRange[1],
                50,0.0,(1.43+0.22*len(self.signals))*max([ymax,1.0])
            )

        axis.GetXaxis().SetLabelSize(0)
        axis.GetXaxis().SetTitle("")
        axis.GetXaxis().SetTickLength(0.015/(1-cv.GetPad(2).GetLeftMargin()-cv.GetPad(2).GetRightMargin()))
        axis.GetYaxis().SetTickLength(0.015/(1-cv.GetPad(2).GetTopMargin()-cv.GetPad(2).GetBottomMargin()))
        #axis.GetYaxis().SetNoExponent(True)
        axis.Draw("AXIS")

        mcStackNominal.Draw("HISTSame")

        rootObj=[]

        for ibin in range(mcHistSumNominal.GetNbinsX()):
            c = mcHistSumNominal.GetBinCenter(ibin+1)
            w = mcHistSumNominal.GetBinWidth(ibin+1)
            m = mcHistSumNominal.GetBinContent(ibin+1)
            err = mcHistSumNominal.GetBinError(ibin+1)

            if m>0.0:
                box = ROOT.TBox(c-0.5*w,m-err,c+0.5*w,m+err)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                box.SetLineWidth(2)
                rootObj.append(box)
                box.Draw("SameF")
                '''
                box2 = ROOT.TBox(c-0.5*w,m-err,c+0.5*w,m+err)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                box2.SetLineWidth(2)
                rootObj.append(box2)
                box2.Draw("SameL")
                '''
                
        for signal in self.signals:
            histSignal = self.getMCSignal(signal)
            histSignal.Draw("HISTSame")
                
        ROOT.gPad.RedrawAxis()
                
        cv.cd(1)
        axisRes=ROOT.TH2F(
            "axisRes"+str(random.randint(0,99999)),";"+xaxisTitle+";Data/Pred.",
            50,self.binRange[0],self.binRange[-1],50,0.1,1.9)#0.55,1.45)#0.1,1.9)
        axisRes.GetYaxis().SetNdivisions(406)
        axisRes.GetXaxis().SetTickLength(0.020/(1-cv.GetPad(1).GetLeftMargin()-cv.GetPad(1).GetRightMargin()))
        axisRes.GetYaxis().SetTickLength(0.015/(1-cv.GetPad(1).GetTopMargin()-cv.GetPad(1).GetBottomMargin()))

        axisRes.Draw("AXIS")

        for ibin in range(mcHistSumNominal.GetNbinsX()):
            c = mcHistSumNominal.GetBinCenter(ibin+1)
            w = mcHistSumNominal.GetBinWidth(ibin+1)
            m = mcHistSumNominal.GetBinContent(ibin+1)

            if m>0.0:
                h = min(mcHistSumNominal.GetBinError(ibin+1)/m,0.899)
                box = ROOT.TBox(c-0.5*w,1-h,c+0.5*w,1+h)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                box.SetLineWidth(2)
                rootObj.append(box)
                box.Draw("SameF")
                box2 = ROOT.TBox(c-0.5*w,1-h,c+0.5*w,1+h)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                box2.SetLineWidth(2)
                rootObj.append(box2)
                box2.Draw("SameL")
        '''
        sumHistRes=sumHistData.Clone(sumHistData.GetName()+str(random.randint(0,99999)))
        rootObj.append(sumHistRes)
        for ibin in range(sumHistRes.GetNbinsX()):
            m = sumHistMC.GetBinContent(ibin+1)
            d = sumHistRes.GetBinContent(ibin+1)
            e = sumHistRes.GetBinError(ibin+1)
            if m>0.0:
                sumHistRes.SetBinContent(ibin+1,d/m)
                sumHistRes.SetBinError(ibin+1,e/m)
            else:
                sumHistRes.SetBinContent(ibin+1,0.0)
                sumHistRes.SetBinError(ibin+1,0)
        sumHistRes.Draw("PESame")
        '''
        axisLine = ROOT.TF1("axisLine"+str(random.randint(0,99999)),"1",self.binRange[0],self.binRange[1])
        axisLine.SetLineColor(ROOT.kBlack)
        axisLine.SetLineWidth(1)
        axisLine.Draw("SameL")
        
        ROOT.gStyle.SetLineWidth(1)
        ROOT.gPad.RedrawAxis()


        hidePave=ROOT.TPaveText(cvxmin-0.06,resHeight-0.029,cvxmin-0.005,resHeight+0.028,"NDC")
        hidePave.SetFillColor(ROOT.kWhite)
        hidePave.SetFillStyle(1001)
        hidePave.Draw("Same")

                
                
        cv.Update()      
        cv.Print(self.plot+self.outputSuffix+".pdf")


bdt = Plot("bdt_SR","BDT discriminant",binRange=[0,1])
bdt.addSignal("HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03_all", 1e5, "", style=[3,2,ROOT.kRed+1])
bdt()
        





'''





f = ROOT.TFile("hists/hist_2016_bdt_SR_nominal.root")

stack = ROOT.THStack()
colors = [ROOT.kOrange,ROOT.kGreen,ROOT.kAzure,ROOT.kMagenta,ROOT.kGray]
for i,sample in enumerate(['topbkg','wjets','dyjets','vgamma','qcd']):
    hist = f.Get(sample)
    print (sample,hist.Integral())
    hist = hist.ProjectionX(hist.GetName()+"tot")
    hist.SetDirectory(0)
    hist.SetFillColor(colors[i])
    #hist.Rebin(2)
    stack.Add(hist,"HISTF")
  
dataSum = None
for i,sample in enumerate(["muon","electron"]):
    hist = f.Get(sample)
    print (sample,hist.Integral())
    hist = hist.ProjectionX(hist.GetName()+"tot")
    hist.SetDirectory(0)
    if dataSum==None:
        dataSum=hist  
    else:
        dataSum.Add(hist)
       
        

signalHist = f.Get('HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03_all')
print ('signal',signalHist.Integral())
signalHist = signalHist.ProjectionX(signalHist.GetName()+"tot")
signalHist.SetDirectory(0)
signalHist.Scale(4e4)
#signalHist.Scale(100)
signalHist.SetLineStyle(2)
signalHist.SetFillStyle(0)
signalHist.SetLineWidth(2)
signalHist.SetLineColor(ROOT.kRed+1)
    
dataSum.SetMarkerStyle(20)
dataSum.SetMarkerSize(1.2)
#dataSum.Rebin(2)

#axis = ROOT.TH2F("axis",";;",50,20,200,50,0.7,1.3*dataSum.GetMaximum())
#axis = ROOT.TH2F("axis",";;",50,40,100,50,0.7,10**(1.2*math.log10(dataSum.GetMaximum())))

axis = ROOT.TH2F("axis",";;",50,0,1,50,0.7,1.3*dataSum.GetMaximum())
#axis = ROOT.TH2F("axis",";;",50,0,1,50,0.7,10**(1.2*math.log10(dataSum.GetMaximum())))
axis.Draw("AXIS") 
stack.Draw("HISTSame")
signalHist.Draw("HISTSame")
dataSum.Draw("PESAME")
#cv.SetLogy(1)
ROOT.gPad.RedrawAxis()
cv.Print("test.pdf")
'''
