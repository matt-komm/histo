import sys
sys.path.append('/vols/cms/vc1117/AN-19-207/classes')
import style
import ROOT
import numpy as np
from array import array

color1 = ROOT.kOrange+7
color2 = ROOT.kSpring+3
color3 = ROOT.kAzure-7

rdf_1 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/15Mar21_2l/2016/HNL_majorana_all_ctau1p0e02_massHNL1p0_Vall5p274e-02-2016/nano_*_Friend.root")
rdf_6 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/15Mar21_2l/2016/HNL_majorana_all_ctau1p0e01_massHNL6p0_Vall1p454e-03-2016/nano_*_Friend.root")
rdf_12 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/15Mar21_2l/2016/HNL_majorana_all_ctau1p0e-01_massHNL12p0_Vall2p314e-03-2016/nano_*_Friend.root")

rdf_1 = rdf_1.Define("subleadingLeptons_dxy_abs", "abs(subleadingLeptons_dxy)")
rdf_6 = rdf_6.Define("subleadingLeptons_dxy_abs", "abs(subleadingLeptons_dxy)")
rdf_12 = rdf_12.Define("subleadingLeptons_dxy_abs", "abs(subleadingLeptons_dxy)")


rdf_emu_coupling_1 = rdf_1.Filter("LHEWeights_coupling_7>0")
rdf_e_coupling_1 = rdf_1.Filter("LHEWeights_coupling_2>0")
rdf_mu_coupling_1 = rdf_1.Filter("LHEWeights_coupling_12>0")
rdf_tau_coupling_1 = rdf_1.Filter("LHEWeights_coupling_67>0")

rdf_emu_coupling_6 = rdf_6.Filter("LHEWeights_coupling_7>0")
rdf_e_coupling_6 = rdf_6.Filter("LHEWeights_coupling_2>0")
rdf_mu_coupling_6 = rdf_6.Filter("LHEWeights_coupling_12>0")
rdf_tau_coupling_6 = rdf_6.Filter("LHEWeights_coupling_67>0")

rdf_emu_coupling_12 = rdf_12.Filter("LHEWeights_coupling_7>0")
rdf_e_coupling_12 = rdf_12.Filter("LHEWeights_coupling_2>0")
rdf_mu_coupling_12 = rdf_12.Filter("LHEWeights_coupling_12>0")
rdf_tau_coupling_12 = rdf_12.Filter("LHEWeights_coupling_67>0")

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_met_emu_coupling_1 = rdf_emu_coupling_1.Histo1D(("hist_met_emu_coupling_1", "hist_met_emu_coupling_1", 50, 0., 150.), "nominal_met")
hist_met_emu_coupling_6 = rdf_emu_coupling_6.Histo1D(("hist_met_emu_coupling_6", "hist_met_emu_coupling_6", 50, 0., 150.), "nominal_met")
hist_met_emu_coupling_12 = rdf_emu_coupling_12.Histo1D(("hist_met_emu_coupling_12", "hist_met_emu_coupling_12", 50, 0., 150.), "nominal_met")

hist_met_tau_coupling_1 = rdf_tau_coupling_1.Histo1D(("hist_met_tau_coupling_1", "hist_met_tau_coupling_1", 25, 0., 150.), "nominal_met")
hist_met_tau_coupling_6 = rdf_tau_coupling_6.Histo1D(("hist_met_tau_coupling_6", "hist_met_tau_coupling_6", 25, 0., 150.), "nominal_met")
hist_met_tau_coupling_12 = rdf_tau_coupling_12.Histo1D(("hist_met_tau_coupling_12", "hist_met_tau_coupling_12", 25, 0., 150.), "nominal_met")

hist_met_emu_coupling_6.SetLineStyle(2)
hist_met_emu_coupling_12.SetLineStyle(3)
hist_met_tau_coupling_6.SetLineStyle(2)
hist_met_tau_coupling_12.SetLineStyle(3)

hist_met_emu_coupling_1.SetLineWidth(2)
hist_met_emu_coupling_6.SetLineWidth(2)
hist_met_emu_coupling_12.SetLineWidth(2)

hist_met_tau_coupling_1.SetLineWidth(2)
hist_met_tau_coupling_6.SetLineWidth(2)
hist_met_tau_coupling_12.SetLineWidth(2)

hist_met_emu_coupling_1.Scale(1./hist_met_emu_coupling_1.Integral())
hist_met_emu_coupling_6.Scale(1./hist_met_emu_coupling_6.Integral())
hist_met_emu_coupling_12.Scale(1./hist_met_emu_coupling_12.Integral())

hist_met_tau_coupling_1.Scale(1./hist_met_tau_coupling_1.Integral())
hist_met_tau_coupling_6.Scale(1./hist_met_tau_coupling_6.Integral())
hist_met_tau_coupling_12.Scale(1./hist_met_tau_coupling_12.Integral())

hist_met_emu_coupling_1.SetLineColor(color1)
hist_met_emu_coupling_6.SetLineColor(color2)
hist_met_emu_coupling_12.SetLineColor(color3)

hist_met_tau_coupling_1.SetLineColor(color1)
hist_met_tau_coupling_6.SetLineColor(color2)
hist_met_tau_coupling_12.SetLineColor(color3)

hist_met_emu_coupling_1.GetXaxis().SetTitle("p^{miss}_{T} [GeV]")
hist_met_emu_coupling_1.GetYaxis().SetTitle("Normalised entries")
hist_met_tau_coupling_1.GetXaxis().SetTitle("p^{miss}_{T} [GeV]")
hist_met_tau_coupling_1.GetYaxis().SetTitle("Normalised entries")


hist_met_emu_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_met_emu_coupling_1.GetMaximum(), hist_met_emu_coupling_6.GetMaximum(), hist_met_emu_coupling_12.GetMaximum()))

hist_met_emu_coupling_1.Draw("HIST")
hist_met_emu_coupling_6.Draw("HIST SAME")
hist_met_emu_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_met_emu_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_met_emu_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_met_emu_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.36, 0.68, 0.77, 0.73, "V_{e} = V_{#mu}, V_{#tau} = 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/met_emu.pdf")

hist_met_tau_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_met_tau_coupling_1.GetMaximum(), hist_met_tau_coupling_6.GetMaximum(), hist_met_tau_coupling_12.GetMaximum()))


hist_met_tau_coupling_1.Draw("HIST")
hist_met_tau_coupling_6.Draw("HIST SAME")
hist_met_tau_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_met_tau_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_met_tau_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_met_tau_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.36, 0.68, 0.77, 0.73, "V_{e} = V_{#mu} = 0, V_{#tau} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")

cv.SaveAs("kinematics/met_tau.pdf")


cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_ht_emu_coupling_1 = rdf_emu_coupling_1.Histo1D(("hist_ht_emu_coupling_1", "hist_ht_emu_coupling_1", 50, 15., 200.), "nominal_ht")
hist_ht_emu_coupling_6 = rdf_emu_coupling_6.Histo1D(("hist_ht_emu_coupling_6", "hist_ht_emu_coupling_6", 50, 15., 200.), "nominal_ht")
hist_ht_emu_coupling_12 = rdf_emu_coupling_12.Histo1D(("hist_ht_emu_coupling_12", "hist_ht_emu_coupling_12", 50, 15., 200.), "nominal_ht")

hist_ht_tau_coupling_1 = rdf_tau_coupling_1.Histo1D(("hist_ht_tau_coupling_1", "hist_ht_tau_coupling_1", 25, 15., 200.), "nominal_ht")
hist_ht_tau_coupling_6 = rdf_tau_coupling_6.Histo1D(("hist_ht_tau_coupling_6", "hist_ht_tau_coupling_6", 25, 15., 200.), "nominal_ht")
hist_ht_tau_coupling_12 = rdf_tau_coupling_12.Histo1D(("hist_ht_tau_coupling_12", "hist_ht_tau_coupling_12", 25, 15., 200.), "nominal_ht")

hist_ht_emu_coupling_6.SetLineStyle(2)
hist_ht_emu_coupling_12.SetLineStyle(3)
hist_ht_tau_coupling_6.SetLineStyle(2)
hist_ht_tau_coupling_12.SetLineStyle(3)

hist_ht_emu_coupling_1.SetLineWidth(2)
hist_ht_emu_coupling_6.SetLineWidth(2)
hist_ht_emu_coupling_12.SetLineWidth(2)

hist_ht_tau_coupling_1.SetLineWidth(2)
hist_ht_tau_coupling_6.SetLineWidth(2)
hist_ht_tau_coupling_12.SetLineWidth(2)

hist_ht_emu_coupling_1.Scale(1./hist_ht_emu_coupling_1.Integral())
hist_ht_emu_coupling_6.Scale(1./hist_ht_emu_coupling_6.Integral())
hist_ht_emu_coupling_12.Scale(1./hist_ht_emu_coupling_12.Integral())

hist_ht_tau_coupling_1.Scale(1./hist_ht_tau_coupling_1.Integral())
hist_ht_tau_coupling_6.Scale(1./hist_ht_tau_coupling_6.Integral())
hist_ht_tau_coupling_12.Scale(1./hist_ht_tau_coupling_12.Integral())

hist_ht_emu_coupling_1.SetLineColor(color1)
hist_ht_emu_coupling_6.SetLineColor(color2)
hist_ht_emu_coupling_12.SetLineColor(color3)

hist_ht_tau_coupling_1.SetLineColor(color1)
hist_ht_tau_coupling_6.SetLineColor(color2)
hist_ht_tau_coupling_12.SetLineColor(color3)

hist_ht_emu_coupling_1.GetXaxis().SetTitle("H_{T} [GeV]")
hist_ht_emu_coupling_1.GetYaxis().SetTitle("Normalised entries")
hist_ht_tau_coupling_1.GetXaxis().SetTitle("H_{T} [GeV]")
hist_ht_tau_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_ht_emu_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_ht_emu_coupling_1.GetMaximum(), hist_ht_emu_coupling_6.GetMaximum(), hist_ht_emu_coupling_12.GetMaximum()))

hist_ht_emu_coupling_1.Draw("HIST")
hist_ht_emu_coupling_6.Draw("HIST SAME")
hist_ht_emu_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_ht_emu_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_ht_emu_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_ht_emu_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.36, 0.68, 0.77, 0.73, "V_{e} = V_{#mu}, V_{#tau} = 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/ht_emu.pdf")

hist_ht_tau_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_ht_tau_coupling_1.GetMaximum(), hist_ht_tau_coupling_6.GetMaximum(), hist_ht_tau_coupling_12.GetMaximum()))


hist_ht_tau_coupling_1.Draw("HIST")
hist_ht_tau_coupling_6.Draw("HIST SAME")
hist_ht_tau_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_ht_tau_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_ht_tau_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_ht_tau_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.36, 0.68, 0.77, 0.73, "V_{e} = V_{#mu} = 0, V_{#tau} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")

cv.SaveAs("kinematics/ht_tau.pdf")



cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_njets_emu_coupling_1 = rdf_emu_coupling_1.Histo1D(("hist_njets_emu_coupling_1", "hist_njets_emu_coupling_1", 5, 0.5, 5.5), "nselectedJets_nominal")
hist_njets_emu_coupling_6 = rdf_emu_coupling_6.Histo1D(("hist_njets_emu_coupling_6", "hist_njets_emu_coupling_6", 5, 0.5, 5.5), "nselectedJets_nominal")
hist_njets_emu_coupling_12 = rdf_emu_coupling_12.Histo1D(("hist_njets_emu_coupling_12", "hist_njets_emu_coupling_12", 5, 0.5, 5.5), "nselectedJets_nominal")

hist_njets_tau_coupling_1 = rdf_tau_coupling_1.Histo1D(("hist_njets_tau_coupling_1", "hist_njets_tau_coupling_1", 5, 0.5, 5.5), "nselectedJets_nominal")
hist_njets_tau_coupling_6 = rdf_tau_coupling_6.Histo1D(("hist_njets_tau_coupling_6", "hist_njets_tau_coupling_6", 5, 0.5, 5.5), "nselectedJets_nominal")
hist_njets_tau_coupling_12 = rdf_tau_coupling_12.Histo1D(("hist_njets_tau_coupling_12", "hist_njets_tau_coupling_12", 5, 0.5, 5.5), "nselectedJets_nominal")

hist_njets_emu_coupling_6.SetLineStyle(2)
hist_njets_emu_coupling_12.SetLineStyle(3)
hist_njets_tau_coupling_6.SetLineStyle(2)
hist_njets_tau_coupling_12.SetLineStyle(3)

hist_njets_emu_coupling_1.SetLineWidth(2)
hist_njets_emu_coupling_6.SetLineWidth(2)
hist_njets_emu_coupling_12.SetLineWidth(2)

hist_njets_tau_coupling_1.SetLineWidth(2)
hist_njets_tau_coupling_6.SetLineWidth(2)
hist_njets_tau_coupling_12.SetLineWidth(2)

hist_njets_emu_coupling_1.Scale(1./hist_njets_emu_coupling_1.Integral())
hist_njets_emu_coupling_6.Scale(1./hist_njets_emu_coupling_6.Integral())
hist_njets_emu_coupling_12.Scale(1./hist_njets_emu_coupling_12.Integral())

hist_njets_tau_coupling_1.Scale(1./hist_njets_tau_coupling_1.Integral())
hist_njets_tau_coupling_6.Scale(1./hist_njets_tau_coupling_6.Integral())
hist_njets_tau_coupling_12.Scale(1./hist_njets_tau_coupling_12.Integral())

hist_njets_emu_coupling_1.SetLineColor(color1)
hist_njets_emu_coupling_6.SetLineColor(color2)
hist_njets_emu_coupling_12.SetLineColor(color3)

hist_njets_tau_coupling_1.SetLineColor(color1)
hist_njets_tau_coupling_6.SetLineColor(color2)
hist_njets_tau_coupling_12.SetLineColor(color3)

hist_njets_emu_coupling_1.GetXaxis().SetTitle("Number of jets")
hist_njets_emu_coupling_1.GetYaxis().SetTitle("Normalised entries")
hist_njets_tau_coupling_1.GetXaxis().SetTitle("Number of jets")
hist_njets_tau_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_njets_emu_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_njets_emu_coupling_1.GetMaximum(), hist_njets_emu_coupling_6.GetMaximum(), hist_njets_emu_coupling_12.GetMaximum()))

hist_njets_emu_coupling_1.Draw("HIST")
hist_njets_emu_coupling_6.Draw("HIST SAME")
hist_njets_emu_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_njets_emu_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_njets_emu_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_njets_emu_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.36, 0.68, 0.77, 0.73, "V_{e} = V_{#mu}, V_{#tau} = 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/njets_emu.pdf")

hist_njets_tau_coupling_1.GetYaxis().SetRangeUser(0, 1.9*max(hist_njets_tau_coupling_1.GetMaximum(), hist_njets_tau_coupling_6.GetMaximum(), hist_njets_tau_coupling_12.GetMaximum()))


hist_njets_tau_coupling_1.Draw("HIST")
hist_njets_tau_coupling_6.Draw("HIST SAME")
hist_njets_tau_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.36, 0.5, 0.77, 0.63)
leg.AddEntry("hist_njets_tau_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_njets_tau_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_njets_tau_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.36, 0.68, 0.77, 0.73, "V_{e} = V_{#mu} = 0, V_{#tau} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")

cv.SaveAs("kinematics/njets_tau.pdf")


cv.Draw("")
hist_dilepton_mass_emu_coupling_1 = rdf_emu_coupling_1.Filter("nsubleadingLeptons==1").Histo1D(("hist_dilepton_mass_emu_coupling_1", "hist_dilepton_mass_emu_coupling_1", 50, 10., 200.), "dilepton_mass")
hist_dilepton_mass_emu_coupling_6 = rdf_emu_coupling_6.Filter("nsubleadingLeptons==1").Histo1D(("hist_dilepton_mass_emu_coupling_6", "hist_dilepton_mass_emu_coupling_6", 50, 10., 200.), "dilepton_mass")
hist_dilepton_mass_emu_coupling_12 = rdf_emu_coupling_12.Filter("nsubleadingLeptons==1").Histo1D(("hist_dilepton_mass_emu_coupling_12", "hist_dilepton_mass_emu_coupling_12", 50, 10., 200.), "dilepton_mass")

hist_dilepton_mass_emu_coupling_6.SetLineStyle(2)
hist_dilepton_mass_emu_coupling_12.SetLineStyle(3)

hist_dilepton_mass_emu_coupling_1.SetLineWidth(2)
hist_dilepton_mass_emu_coupling_6.SetLineWidth(2)
hist_dilepton_mass_emu_coupling_12.SetLineWidth(2)

hist_dilepton_mass_emu_coupling_1.Scale(1./hist_dilepton_mass_emu_coupling_1.Integral())
hist_dilepton_mass_emu_coupling_6.Scale(1./hist_dilepton_mass_emu_coupling_6.Integral())
hist_dilepton_mass_emu_coupling_12.Scale(1./hist_dilepton_mass_emu_coupling_12.Integral())

hist_dilepton_mass_emu_coupling_1.SetLineColor(color1)
hist_dilepton_mass_emu_coupling_6.SetLineColor(color2)
hist_dilepton_mass_emu_coupling_12.SetLineColor(color3)


hist_dilepton_mass_emu_coupling_1.GetXaxis().SetTitle("m(l_{1}, l_{2}) [GeV]")
hist_dilepton_mass_emu_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_dilepton_mass_emu_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_dilepton_mass_emu_coupling_1.GetMaximum(), hist_dilepton_mass_emu_coupling_6.GetMaximum(), hist_dilepton_mass_emu_coupling_12.GetMaximum()))

hist_dilepton_mass_emu_coupling_1.Draw("HIST")
hist_dilepton_mass_emu_coupling_6.Draw("HIST SAME")
hist_dilepton_mass_emu_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_dilepton_mass_emu_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_dilepton_mass_emu_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_dilepton_mass_emu_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{e} = V_{#mu}, V_{#tau} = 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/dilepton_mass_emu.pdf")

cv.Draw("")
hist_dRl2j_coupling_1 = rdf_emu_coupling_1.Filter("nsubleadingLeptons==1").Histo1D(("hist_dRl2j_coupling_1", "hist_dRl2j_coupling_1", 50, 0., 3.), "nominal_dR_l2j")
hist_dRl2j_coupling_6 = rdf_emu_coupling_6.Filter("nsubleadingLeptons==1").Histo1D(("hist_dRl2j_coupling_6", "hist_dRl2j_coupling_6", 50, 0., 3.), "nominal_dR_l2j")
hist_dRl2j_coupling_12 = rdf_emu_coupling_12.Filter("nsubleadingLeptons==1").Histo1D(("hist_dRl2j_coupling_12", "hist_dRl2j_coupling_12", 50, 0., 3.), "nominal_dR_l2j")

hist_dRl2j_coupling_6.SetLineStyle(2)
hist_dRl2j_coupling_12.SetLineStyle(3)

hist_dRl2j_coupling_1.SetLineWidth(2)
hist_dRl2j_coupling_6.SetLineWidth(2)
hist_dRl2j_coupling_12.SetLineWidth(2)

hist_dRl2j_coupling_1.Scale(1./hist_dRl2j_coupling_1.Integral())
hist_dRl2j_coupling_6.Scale(1./hist_dRl2j_coupling_6.Integral())
hist_dRl2j_coupling_12.Scale(1./hist_dRl2j_coupling_12.Integral())

hist_dRl2j_coupling_1.SetLineColor(color1)
hist_dRl2j_coupling_6.SetLineColor(color2)
hist_dRl2j_coupling_12.SetLineColor(color3)


hist_dRl2j_coupling_1.GetXaxis().SetTitle("min#DeltaR(l_{2}, jets)")
hist_dRl2j_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_dRl2j_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_dRl2j_coupling_1.GetMaximum(), hist_dRl2j_coupling_6.GetMaximum(), hist_dRl2j_coupling_12.GetMaximum()))

hist_dRl2j_coupling_1.Draw("HIST")
hist_dRl2j_coupling_6.Draw("HIST SAME")
hist_dRl2j_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_dRl2j_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_dRl2j_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_dRl2j_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{e} = V_{#mu}, V_{#tau} = 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/dR_l2j_emu.pdf")



hist_dilepton_mass_tau_coupling_1 = rdf_tau_coupling_1.Filter("nsubleadingLeptons==1").Histo1D(("hist_dilepton_mass_tau_coupling_1", "hist_dilepton_mass_tau_coupling_1", 50, 10., 200.), "dilepton_mass")
hist_dilepton_mass_tau_coupling_6 = rdf_tau_coupling_6.Filter("nsubleadingLeptons==1").Histo1D(("hist_dilepton_mass_tau_coupling_6", "hist_dilepton_mass_tau_coupling_6", 50, 10., 200.), "dilepton_mass")
hist_dilepton_mass_tau_coupling_12 = rdf_tau_coupling_12.Filter("nsubleadingLeptons==1").Histo1D(("hist_dilepton_mass_tau_coupling_12", "hist_dilepton_mass_tau_coupling_12", 50, 10., 200.), "dilepton_mass")

hist_dilepton_mass_tau_coupling_6.SetLineStyle(2)
hist_dilepton_mass_tau_coupling_12.SetLineStyle(3)

hist_dilepton_mass_tau_coupling_1.SetLineWidth(2)
hist_dilepton_mass_tau_coupling_6.SetLineWidth(2)
hist_dilepton_mass_tau_coupling_12.SetLineWidth(2)

hist_dilepton_mass_tau_coupling_1.Scale(1./hist_dilepton_mass_tau_coupling_1.Integral())
hist_dilepton_mass_tau_coupling_6.Scale(1./hist_dilepton_mass_tau_coupling_6.Integral())
hist_dilepton_mass_tau_coupling_12.Scale(1./hist_dilepton_mass_tau_coupling_12.Integral())

hist_dilepton_mass_tau_coupling_1.SetLineColor(color1)
hist_dilepton_mass_tau_coupling_6.SetLineColor(color2)
hist_dilepton_mass_tau_coupling_12.SetLineColor(color3)


hist_dilepton_mass_tau_coupling_1.GetXaxis().SetTitle("m(l_{1}, l_{2}) [GeV]")
hist_dilepton_mass_tau_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_dilepton_mass_tau_coupling_1.GetYaxis().SetRangeUser(0, 1.4*max(hist_dilepton_mass_tau_coupling_1.GetMaximum(), hist_dilepton_mass_tau_coupling_6.GetMaximum(), hist_dilepton_mass_tau_coupling_12.GetMaximum()))

hist_dilepton_mass_tau_coupling_1.Draw("HIST")
hist_dilepton_mass_tau_coupling_6.Draw("HIST SAME")
hist_dilepton_mass_tau_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_dilepton_mass_tau_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_dilepton_mass_tau_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_dilepton_mass_tau_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{e} = V_{#mu} = 0, V_{#tau} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/dilepton_mass_tau.pdf")

rdf_e_coupling_1_electron = rdf_e_coupling_1.Filter("nsubleadingLeptons==1").Filter("subleadingLeptons_isElectron[0] == 1")
rdf_e_coupling_6_electron = rdf_e_coupling_6.Filter("nsubleadingLeptons==1").Filter("subleadingLeptons_isElectron[0] == 1")
rdf_e_coupling_12_electron = rdf_e_coupling_12.Filter("nsubleadingLeptons==1").Filter("subleadingLeptons_isElectron[0] == 1")

rdf_mu_coupling_1_muon = rdf_mu_coupling_1.Filter("nsubleadingLeptons==1").Filter("subleadingLeptons_isMuon[0] == 1")
rdf_mu_coupling_6_muon = rdf_mu_coupling_6.Filter("nsubleadingLeptons==1").Filter("subleadingLeptons_isMuon[0] == 1")
rdf_mu_coupling_12_muon = rdf_mu_coupling_12.Filter("nsubleadingLeptons==1").Filter("subleadingLeptons_isMuon[0] == 1")

cv.Draw("")
cv.SetLogx()
nbins_dxy = 25
bins_dxy = array('d', np.logspace(-3., 2., num=nbins_dxy+1))

hist_subleadingLeptons_dxy_e_coupling_1 = rdf_e_coupling_1_electron.Histo1D(("hist_subleadingLeptons_dxy_e_coupling_1", "hist_subleadingLeptons_dxy_e_coupling_1", nbins_dxy, bins_dxy), "subleadingLeptons_dxy_abs")
hist_subleadingLeptons_dxy_e_coupling_6 = rdf_e_coupling_6_electron.Histo1D(("hist_subleadingLeptons_dxy_e_coupling_6", "hist_subleadingLeptons_dxy_e_coupling_6", nbins_dxy, bins_dxy), "subleadingLeptons_dxy_abs")
hist_subleadingLeptons_dxy_e_coupling_12 = rdf_e_coupling_12_electron.Histo1D(("hist_subleadingLeptons_dxy_e_coupling_12", "hist_subleadingLeptons_dxy_e_coupling_12", nbins_dxy, bins_dxy), "subleadingLeptons_dxy_abs")

hist_subleadingLeptons_dxy_e_coupling_6.SetLineStyle(2)
hist_subleadingLeptons_dxy_e_coupling_12.SetLineStyle(3)

hist_subleadingLeptons_dxy_e_coupling_1.SetLineWidth(2)
hist_subleadingLeptons_dxy_e_coupling_6.SetLineWidth(2)
hist_subleadingLeptons_dxy_e_coupling_12.SetLineWidth(2)

hist_subleadingLeptons_dxy_e_coupling_1.Scale(1./hist_subleadingLeptons_dxy_e_coupling_1.Integral())
hist_subleadingLeptons_dxy_e_coupling_6.Scale(1./hist_subleadingLeptons_dxy_e_coupling_6.Integral())
hist_subleadingLeptons_dxy_e_coupling_12.Scale(1./hist_subleadingLeptons_dxy_e_coupling_12.Integral())

hist_subleadingLeptons_dxy_e_coupling_1.SetLineColor(color1)
hist_subleadingLeptons_dxy_e_coupling_6.SetLineColor(color2)
hist_subleadingLeptons_dxy_e_coupling_12.SetLineColor(color3)


hist_subleadingLeptons_dxy_e_coupling_1.GetXaxis().SetTitle("e_{2} |d_{xy}| [cm]")
hist_subleadingLeptons_dxy_e_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_subleadingLeptons_dxy_e_coupling_1.GetYaxis().SetRangeUser(0, 1.7*max(hist_subleadingLeptons_dxy_e_coupling_1.GetMaximum(), hist_subleadingLeptons_dxy_e_coupling_6.GetMaximum(), hist_subleadingLeptons_dxy_e_coupling_12.GetMaximum()))

hist_subleadingLeptons_dxy_e_coupling_1.Draw("HIST")
hist_subleadingLeptons_dxy_e_coupling_6.Draw("HIST SAME")
hist_subleadingLeptons_dxy_e_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_subleadingLeptons_dxy_e_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_subleadingLeptons_dxy_e_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_subleadingLeptons_dxy_e_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{#mu} = V_{#tau} = 0, V_{e} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/subleadingLeptons_dxy_e.pdf")

cv.Draw("")
cv.SetLogx()
nbins_dxy = 25
bins_dxy = array('d', np.logspace(-3., 2., num=nbins_dxy+1))

hist_subleadingLeptons_dxy_mu_coupling_1 = rdf_mu_coupling_1_muon.Histo1D(("hist_subleadingLeptons_dxy_mu_coupling_1", "hist_subleadingLeptons_dxy_mu_coupling_1", nbins_dxy, bins_dxy), "subleadingLeptons_dxy_abs")
hist_subleadingLeptons_dxy_mu_coupling_6 = rdf_mu_coupling_6_muon.Histo1D(("hist_subleadingLeptons_dxy_mu_coupling_6", "hist_subleadingLeptons_dxy_mu_coupling_6", nbins_dxy, bins_dxy), "subleadingLeptons_dxy_abs")
hist_subleadingLeptons_dxy_mu_coupling_12 = rdf_mu_coupling_12_muon.Histo1D(("hist_subleadingLeptons_dxy_mu_coupling_12", "hist_subleadingLeptons_dxy_mu_coupling_12", nbins_dxy, bins_dxy), "subleadingLeptons_dxy_abs")

hist_subleadingLeptons_dxy_mu_coupling_6.SetLineStyle(2)
hist_subleadingLeptons_dxy_mu_coupling_12.SetLineStyle(3)

hist_subleadingLeptons_dxy_mu_coupling_1.SetLineWidth(2)
hist_subleadingLeptons_dxy_mu_coupling_6.SetLineWidth(2)
hist_subleadingLeptons_dxy_mu_coupling_12.SetLineWidth(2)

hist_subleadingLeptons_dxy_mu_coupling_1.Scale(1./hist_subleadingLeptons_dxy_mu_coupling_1.Integral())
hist_subleadingLeptons_dxy_mu_coupling_6.Scale(1./hist_subleadingLeptons_dxy_mu_coupling_6.Integral())
hist_subleadingLeptons_dxy_mu_coupling_12.Scale(1./hist_subleadingLeptons_dxy_mu_coupling_12.Integral())

hist_subleadingLeptons_dxy_mu_coupling_1.SetLineColor(color1)
hist_subleadingLeptons_dxy_mu_coupling_6.SetLineColor(color2)
hist_subleadingLeptons_dxy_mu_coupling_12.SetLineColor(color3)


hist_subleadingLeptons_dxy_mu_coupling_1.GetXaxis().SetTitle("#mu_{2} |d_{xy}| [cm]")
hist_subleadingLeptons_dxy_mu_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_subleadingLeptons_dxy_mu_coupling_1.GetYaxis().SetRangeUser(0, 1.7*max(hist_subleadingLeptons_dxy_mu_coupling_1.GetMaximum(), hist_subleadingLeptons_dxy_mu_coupling_6.GetMaximum(), hist_subleadingLeptons_dxy_mu_coupling_12.GetMaximum()))

hist_subleadingLeptons_dxy_mu_coupling_1.Draw("HIST")
hist_subleadingLeptons_dxy_mu_coupling_6.Draw("HIST SAME")
hist_subleadingLeptons_dxy_mu_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_subleadingLeptons_dxy_mu_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_subleadingLeptons_dxy_mu_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_subleadingLeptons_dxy_mu_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{e} = V_{#tau} = 0, V_{#mu} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/subleadingLeptons_dxy_mu.pdf")


cv.Draw("")
cv.SetLogx(0)

hist_subleadingLeptons_pt_mu_coupling_1 = rdf_mu_coupling_1_muon.Histo1D(("hist_subleadingLeptons_pt_mu_coupling_1", "hist_subleadingLeptons_pt_mu_coupling_1", 20, 5., 40.), "subleadingLeptons_pt")
hist_subleadingLeptons_pt_mu_coupling_6 = rdf_mu_coupling_6_muon.Histo1D(("hist_subleadingLeptons_pt_mu_coupling_6", "hist_subleadingLeptons_pt_mu_coupling_6",  20, 5., 40.), "subleadingLeptons_pt")
hist_subleadingLeptons_pt_mu_coupling_12 = rdf_mu_coupling_12_muon.Histo1D(("hist_subleadingLeptons_pt_mu_coupling_12", "hist_subleadingLeptons_pt_mu_coupling_12",  20, 5., 40.), "subleadingLeptons_pt")

hist_subleadingLeptons_pt_mu_coupling_6.SetLineStyle(2)
hist_subleadingLeptons_pt_mu_coupling_12.SetLineStyle(3)

hist_subleadingLeptons_pt_mu_coupling_1.SetLineWidth(2)
hist_subleadingLeptons_pt_mu_coupling_6.SetLineWidth(2)
hist_subleadingLeptons_pt_mu_coupling_12.SetLineWidth(2)

hist_subleadingLeptons_pt_mu_coupling_1.Scale(1./hist_subleadingLeptons_pt_mu_coupling_1.Integral())
hist_subleadingLeptons_pt_mu_coupling_6.Scale(1./hist_subleadingLeptons_pt_mu_coupling_6.Integral())
hist_subleadingLeptons_pt_mu_coupling_12.Scale(1./hist_subleadingLeptons_pt_mu_coupling_12.Integral())

hist_subleadingLeptons_pt_mu_coupling_1.SetLineColor(color1)
hist_subleadingLeptons_pt_mu_coupling_6.SetLineColor(color2)
hist_subleadingLeptons_pt_mu_coupling_12.SetLineColor(color3)


hist_subleadingLeptons_pt_mu_coupling_1.GetXaxis().SetTitle("#mu_{2} p_{T} [GeV]")
hist_subleadingLeptons_pt_mu_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_subleadingLeptons_pt_mu_coupling_1.GetYaxis().SetRangeUser(0, 1.7*max(hist_subleadingLeptons_pt_mu_coupling_1.GetMaximum(), hist_subleadingLeptons_pt_mu_coupling_6.GetMaximum(), hist_subleadingLeptons_pt_mu_coupling_12.GetMaximum()))

hist_subleadingLeptons_pt_mu_coupling_1.Draw("HIST")
hist_subleadingLeptons_pt_mu_coupling_6.Draw("HIST SAME")
hist_subleadingLeptons_pt_mu_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_subleadingLeptons_pt_mu_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_subleadingLeptons_pt_mu_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_subleadingLeptons_pt_mu_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{e} = V_{#tau} = 0, V_{#mu} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/subleadingLeptons_pt_mu.pdf")

cv.Draw("")

hist_subleadingLeptons_pt_e_coupling_1 = rdf_e_coupling_1_electron.Histo1D(("hist_subleadingLeptons_pt_e_coupling_1", "hist_subleadingLeptons_pt_e_coupling_1", 20, 5., 40.), "subleadingLeptons_pt")
hist_subleadingLeptons_pt_e_coupling_6 = rdf_e_coupling_6_electron.Histo1D(("hist_subleadingLeptons_pt_e_coupling_6", "hist_subleadingLeptons_pt_e_coupling_6",  20, 5., 40.), "subleadingLeptons_pt")
hist_subleadingLeptons_pt_e_coupling_12 = rdf_e_coupling_12_electron.Histo1D(("hist_subleadingLeptons_pt_e_coupling_12", "hist_subleadingLeptons_pt_e_coupling_12",  20, 5., 40.), "subleadingLeptons_pt")

hist_subleadingLeptons_pt_e_coupling_6.SetLineStyle(2)
hist_subleadingLeptons_pt_e_coupling_12.SetLineStyle(3)

hist_subleadingLeptons_pt_e_coupling_1.SetLineWidth(2)
hist_subleadingLeptons_pt_e_coupling_6.SetLineWidth(2)
hist_subleadingLeptons_pt_e_coupling_12.SetLineWidth(2)

hist_subleadingLeptons_pt_e_coupling_1.Scale(1./hist_subleadingLeptons_pt_e_coupling_1.Integral())
hist_subleadingLeptons_pt_e_coupling_6.Scale(1./hist_subleadingLeptons_pt_e_coupling_6.Integral())
hist_subleadingLeptons_pt_e_coupling_12.Scale(1./hist_subleadingLeptons_pt_e_coupling_12.Integral())

hist_subleadingLeptons_pt_e_coupling_1.SetLineColor(color1)
hist_subleadingLeptons_pt_e_coupling_6.SetLineColor(color2)
hist_subleadingLeptons_pt_e_coupling_12.SetLineColor(color3)


hist_subleadingLeptons_pt_e_coupling_1.GetXaxis().SetTitle("e_{2} p_{T} [GeV]")
hist_subleadingLeptons_pt_e_coupling_1.GetYaxis().SetTitle("Normalised entries")

hist_subleadingLeptons_pt_e_coupling_1.GetYaxis().SetRangeUser(0, 1.7*max(hist_subleadingLeptons_pt_e_coupling_1.GetMaximum(), hist_subleadingLeptons_pt_e_coupling_6.GetMaximum(), hist_subleadingLeptons_pt_e_coupling_12.GetMaximum()))

hist_subleadingLeptons_pt_e_coupling_1.Draw("HIST")
hist_subleadingLeptons_pt_e_coupling_6.Draw("HIST SAME")
hist_subleadingLeptons_pt_e_coupling_12.Draw("HIST SAME")

leg = style.makeLegend(0.39, 0.5, 0.80, 0.63)
leg.AddEntry("hist_subleadingLeptons_pt_e_coupling_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_subleadingLeptons_pt_e_coupling_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_subleadingLeptons_pt_e_coupling_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.39, 0.68, 0.80, 0.73, "V_{#mu} = V_{#tau} = 0, V_{e} #neq 0")

style.makeText(0.27, 0.75, 0.4, 0.75, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeCMSText(0.17, 0.85, additionalText="Simulation Preliminary")
cv.SaveAs("kinematics/subleadingLeptons_pt_e.pdf")





