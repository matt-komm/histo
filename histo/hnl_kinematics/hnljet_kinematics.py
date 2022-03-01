from histo import style
import ROOT
import numpy as np
from array import array

color1 = ROOT.kOrange+7
color2 = ROOT.kSpring+3
color3 = ROOT.kAzure-7

rdf = ROOT.RDataFrame("Friends", "/home/hep/vc1117/LLP/CMSSW_10_2_18/src/HNL_dirac_pt20_ctau1p0e-03_massHNL16p0_Vall1p551e-02_*_Friend.root")

rdf = rdf.Filter("nelectron+nmuon>0")
rdf = rdf.Define("selectedJets_nominal_isLLP_QL", "selectedJets_nominal_isLLP_QMU+selectedJets_nominal_isLLP_QE")

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_Q = rdf.Histo1D(("hist_Q", "hist_Q", 50, 15., 60.), "selectedJets_nominal_pt", "selectedJets_nominal_isLLP_Q")
hist_QL = rdf.Histo1D(("hist_QL", "hist_QL", 50, 15., 60.), "selectedJets_nominal_pt", "selectedJets_nominal_isLLP_QL")

hist_QL.SetLineStyle(2)

hist_QL.SetLineWidth(2)


hist_Q.Scale(1./hist_Q.Integral())
hist_QL.Scale(1./hist_QL.Integral())

hist_Q.SetLineColor(color1)
hist_QL.SetLineColor(color2)

hist_Q.GetXaxis().SetTitle("HNL jet p_{T} [GeV]")
hist_Q.GetYaxis().SetTitle("Normalised entries")


hist_Q.GetYaxis().SetRangeUser(0, 1.4*max(hist_Q.GetMaximum(), hist_QL.GetMaximum()))

hist_Q.Draw("HIST")
hist_QL.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_Q", "jet (resolved)", "l")
leg.AddEntry("hist_QL", "jet w/ lepton (merged)", "l")
leg.Draw("")


style.makeText(0.27, 0.75, 0.4, 0.8, "m_{N} = 16 GeV, V_{e} = V_{#mu} = V_{#tau}")
style.makeText(0.27, 0.81, 0.4, 0.86, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowl_{1}N#rightarrowl_{1}l_{2}q#bar{q'}")
#style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("hnljet_pt.pdf")


cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_Q = rdf.Histo1D(("hist_Q", "hist_Q", 50, 0, 2), "selectedJets_nominal_minDeltaRSubtraction", "selectedJets_nominal_isLLP_Q")
hist_QL = rdf.Histo1D(("hist_QL", "hist_QL", 50, 0, 2), "selectedJets_nominal_minDeltaRSubtraction", "selectedJets_nominal_isLLP_QMU")

hist_QL.SetLineStyle(2)

hist_QL.SetLineWidth(2)


hist_Q.Scale(1./hist_Q.Integral())
hist_QL.Scale(1./hist_QL.Integral())

hist_Q.SetLineColor(color1)
hist_QL.SetLineColor(color2)

hist_Q.GetXaxis().SetTitle("#DeltaR(HNL jet, l_{2}) [GeV]")
hist_Q.GetYaxis().SetTitle("Normalised entries")


hist_Q.GetYaxis().SetRangeUser(0, 1.4*max(hist_Q.GetMaximum(), hist_QL.GetMaximum()))

hist_Q.Draw("HIST")
hist_QL.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_Q", "jet (resolved)", "l")
leg.AddEntry("hist_QL", "jet w/ lepton (merged)", "l")
leg.Draw("")

style.makeText(0.27, 0.75, 0.4, 0.8, "m_{N} = 16 GeV, V_{e} = V_{#mu} = V_{#tau}")
style.makeText(0.27, 0.81, 0.4, 0.86, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowl_{1}N#rightarrowl_{1}l_{2}q#bar{q'}")
#style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("hnljet_dR.pdf")


