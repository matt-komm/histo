import ROOT
import numpy as np
import os
from array import array
from histo import style
from histo.lumi import lumi


class Variable:
    """
    This class is responsible for making the histogram and plotting it for a given variable
    """
    def __init__(self, varexp, name, nbins, xmin, xmax, logx=False, logy=True, corrections=False):
        print(varexp, name, nbins, xmin, xmax, logx, logy)
        self.varexp = varexp
        self.stack = ROOT.THStack(varexp, varexp)
        self.corrections = corrections
        if corrections:
            self.stack_up = ROOT.THStack(varexp, varexp)
        self.signals = []
        self.data = None
        self.name = name
        self.logx = logx
        self.logy = logy
        self.leg = style.makeLegend(0.82, 0.50, 0.99, 0.9)
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.leg.SetTextSize(self.leg.GetTextSize()*0.8)
        if self.logx:
            bins = np.geomspace(float(xmin), float(xmax), num=int(nbins)+1)
            bins = array('d', bins)
            self.args = (varexp, varexp, nbins, bins)
        else:
            self.args = (varexp, varexp, nbins, xmin, xmax)

    def Add(self, hist, title, isSignal=False, isData=False, correction="nominal"):
        hist.SetDirectory(0)
        if isSignal:
            hist.SetLineStyle(1+len(self.signals))
            hist.SetTitle(title)
            self.signals.append(hist)
            #self.leg.AddEntry(hist, title, "l")
        elif isData:
            self.data = hist
            self.leg.AddEntry(hist, title, "p")
        else:
            if correction == "nominal":
                self.stack.Add(hist)
                self.leg.AddEntry(hist, title, "f")
            elif correction == "up":
                self.stack_up.Add(hist)
            elif correction == "down":
                self.stack_down.Add(hist)

    def Draw(self, suffix, opt, draw_text, year="2016", output_dir="test"):
        print ("plotting "+self.varexp)
        canvas = style.makeCanvas(name=self.varexp)


        upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
        lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)
        upperPad.SetBottomMargin(0.00001)
        upperPad.SetBorderMode(0)
        upperPad.SetTopMargin(0.11)
        upperPad.SetRightMargin(0.2)
        upperPad.SetLeftMargin(0.13)
        lowerPad.SetTopMargin(0.00001)
        lowerPad.SetBottomMargin(0.4)
        lowerPad.SetRightMargin(0.2)
        lowerPad.SetLeftMargin(0.13)
        lowerPad.SetBorderMode(0)
        #canvas.SetBottomMargin(0.2)
        #canvas.SetTopMargin(0.05)
        upperPad.Draw()
        lowerPad.Draw()
        upperPad.cd()

        self.stack.Draw(opt)
        self.stack.GetYaxis().SetTitle("Entries/bin")
        self.stack.GetYaxis().SetTitleOffset(1.2)
        self.stack.SetMinimum(10)

        for signal in self.signals:
            signal.SetFillStyle(0)
            signal.Draw("HIST SAME")

        if self.logy:
            if len(self.signals) > 0:
                self.stack.SetMaximum(self.stack.GetMaximum()*70)
            else:
                self.stack.SetMaximum(self.stack.GetMaximum()*10)
            self.stack.SetMinimum(10)
            upperPad.SetLogy()
        else:
            self.stack.SetMaximum(self.stack.GetMaximum()*1.5)

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
        axis.GetYaxis().SetTitle("Data/MC")
        axis.GetYaxis().SetTitleOffset(1.2)
        axis.GetXaxis().SetTitleOffset(3.5)
        if self.logx:
            lowerPad.SetLogx()
        axis.Draw("AXIS")


        if self.corrections:
            self.corr_up = self.stack_up.GetStack().Last()


        rootObj = []
        rootObj.append(axis)

        line = ROOT.TLine(self.xmin, 1, self.xmax, 1)
        line.Draw("SAME")
        rootObj.append(line)


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
        self.leg.AddEntry(box, "MC unc.", "f")


        if self.corrections:
            for ibin in range(self.sumMC.GetNbinsX()):
                c = self.sumMC.GetBinCenter(ibin+1)
                w = self.sumMC.GetBinWidth(ibin+1)
                m = self.sumMC.GetBinContent(ibin+1)

                if m > 0.0:
                    h = (self.corr_up.GetBinContent(ibin+1)-m)/m
                    box = ROOT.TBox(c-0.25*w, 1-h, c+0.25*w, 1+h)
                    #box.SetFillStyle(3345)
                    box.SetFillColor(ROOT.TColor.GetColor("#d4caf0"))
                    box.SetLineColor(ROOT.TColor.GetColor("#7170bf"))
                    rootObj.append(box)
                    box.Draw("SameF")
                    box2 = ROOT.TBox(c-0.25*w, 1-h, c+0.25*w, 1+h)
                    box2.SetFillStyle(0)
                    box2.SetFillColor(ROOT.TColor.GetColor("#d4caf0"))
                    box2.SetLineColor(ROOT.TColor.GetColor("#7170bf"))
                    rootObj.append(box2)
                    box2.Draw("SameL")
            self.leg.AddEntry(box, "track unc.", "f")


        if self.data is not None:
            style.makeText(0.13, 0.13, 0.2, 0.2, "Data/MC = {0:.3g}".format(self.data.Integral()/sumMC.Integral()))
            for ibin in range(self.hist_ratio.GetNbinsX()):
                e = self.data.GetBinError(ibin+1)
                m = self.data.GetBinContent(ibin+1)
                if m > 0.0:
                    self.hist_ratio.SetBinError(ibin+1, e/m)
                else:
                    self.hist_ratio.SetBinError(ibin+1, 0)

            self.hist_ratio.Draw("PE SAME")

        canvas.cd()
        self.leg.Draw("SAME")
        if len(self.signals) > 0:
            self.leg_signal = style.makeLegend(0.13, 0.78, 0.67, 0.89)
            self.leg_signal.SetTextSize(25)
            for signal in self.signals:
                self.leg_signal.AddEntry(signal, signal.GetTitle(), "l")
            self.leg_signal.Draw("SAME")
            style.makeLumiText(0.58, 0.91, lumi[year], year, size=25)

        else:
            style.makeLumiText(0.57, 0.91, lumi[year], year, size=25)

        style.makeText(0.13, 0.99, 0.5, 0.99, draw_text, size=25)
        canvas.SaveAs(os.path.join(output_dir, suffix+self.varexp.replace("/","_")+"_"+year+"_thesis.pdf"))

        style.makeCMSText(0.15, 0.91, additionalText="Preliminary")
        canvas.SaveAs(os.path.join(output_dir, suffix+self.varexp.replace("/","_")+"_"+year+".root"))
        canvas.SaveAs(os.path.join(output_dir, suffix+self.varexp.replace("/","_")+"_"+year+".pdf"))
