import ROOT
from . import style
import numpy as np
import os
from array import array

lumi = {"2016": 35.92, "2017": 41.53, "2018": 59.68}


# This class is responsible for making the histogram and plotting it for a given variable
class Variable:
    def __init__(self, varexp, name, nbins, xmin, xmax, logx=False, logy=True):
        print(varexp, name, nbins, xmin, xmax, logx, logy)
        self.varexp = varexp
        self.stack = ROOT.THStack(varexp, varexp)
        self.signals = []
        self.data = None
        self.name = name
        self.logx = logx
        self.logy = logy
        self.leg = style.makeLegend(0.65, 0.70, 0.90, 0.88)
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.leg.SetTextSize(self.leg.GetTextSize()*0.8)
        if self.logx:
            bins = np.geomspace(float(xmin), float(xmax), num=int(nbins)+1)
            bins = array('d', bins)
            self.args = (varexp, varexp, nbins, bins)
        else:
            self.args = (varexp, varexp, nbins, xmin, xmax)

    def Add(self, hist, title, isSignal=False, isData=False):
        hist.SetDirectory(0)
        if isSignal:
            self.signals.append(hist)
            hist.SetLineStyle(1+len(self.signals))
            self.leg.AddEntry(hist, title, "l")
        elif isData:
            self.data = hist
            self.leg.AddEntry(hist, title, "p")
        else:
            self.stack.Add(hist)
            self.leg.AddEntry(hist, title, "f")
    def Draw(self, suffix, opt, draw_text, year="2016", output_dir="test"):
        print ("plotting "+self.varexp)
        canvas = style.makeCanvas(name=self.varexp)


        upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
        lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)
        upperPad.SetBottomMargin(0.00001)
        upperPad.SetBorderMode(0)
        upperPad.SetTopMargin(0.15)
        lowerPad.SetTopMargin(0.00001)
        lowerPad.SetBottomMargin(0.4)
        lowerPad.SetBorderMode(0)
        canvas.SetBottomMargin(0.2)
        canvas.SetTopMargin(0.1)
        upperPad.Draw()
        lowerPad.Draw()
        upperPad.cd()

        self.stack.Draw(opt)
        self.stack.SetMinimum(1)

        for signal in self.signals:
            signal.SetFillStyle(0)
            signal.Draw("HIST SAME")

        if self.logy:
            self.stack.SetMaximum(self.stack.GetMaximum()*100)
            upperPad.SetLogy()
        else:
            self.stack.SetMaximum(self.stack.GetMaximum()*2)

        if self.logx:
            upperPad.SetLogx()

        sumMC = self.stack.GetStack().Last()
        self.sumMC = sumMC.Clone("sum")

        if self.data is not None:
            self.data.Draw("P SAME")
            self.hist_ratio = self.data.Clone("ratio histogram")
            self.hist_ratio.Divide(self.sumMC)

        lowerPad.cd()
        axis = self.sumMC.Clone("axis")
        axis.SetMinimum(0.25)
        axis.SetMaximum(1.75)
        axis.GetXaxis().SetTitle(self.name)
        axis.GetXaxis().SetTitleOffset(2.5)
        if self.logx:
            lowerPad.SetLogx()
        axis.Draw("AXIS")

        rootObj = []
        rootObj.append(axis)

        line = ROOT.TLine(self.xmin, 1, self.xmax, 1)
        line.Draw("SAME")
        rootObj.append(line)

        if self.data is not None:
            style.makeText(0.1, 0.2, 0.2, 0.3, "data/MC = {0:.3g}".format(self.data.Integral()/sumMC.Integral()))
            for ibin in range(self.hist_ratio.GetNbinsX()):
                e = self.data.GetBinError(ibin+1)
                m = self.data.GetBinContent(ibin+1)
                if m > 0.0:
                    self.hist_ratio.SetBinError(ibin+1, e/m)
                else:
                    self.hist_ratio.SetBinError(ibin+1, 0)

            self.hist_ratio.Draw("PE SAME")

        for ibin in range(self.sumMC.GetNbinsX()):
            c = self.sumMC.GetBinCenter(ibin+1)
            w = self.sumMC.GetBinWidth(ibin+1)
            m = self.sumMC.GetBinContent(ibin+1)
            if m > 0.0:
                h = min(self.sumMC.GetBinError(ibin+1)/m, 0.399)
                box = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                rootObj.append(box)
                box.Draw("SameF")
                box2 = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                rootObj.append(box2)
                box2.Draw("SameL")

        canvas.cd()
        self.leg.Draw("SAME")
        style.makeCMSText(0.13, 0.88, additionalText="Preliminary")
        style.makeText(0.15, 0.94, 0.5, 0.94, draw_text)
        style.makeLumiText(0.8, 0.97, lumi[year], year)
        canvas.SaveAs(os.path.join(output_dir, suffix+self.varexp.replace("/","_")+"_"+year+".pdf"))
        canvas.SaveAs(os.path.join(output_dir, suffix+self.varexp.replace("/","_")+"_"+year+".png"))
