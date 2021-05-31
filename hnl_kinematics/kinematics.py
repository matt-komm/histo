from histo import style
import ROOT
import numpy as np

color1 = ROOT.kOrange+7
color2 = ROOT.kSpring+3
color3 = ROOT.kAzure-7

rdf = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/leptonEff/03Nov20/HNL_majorana_all_ctau1p0e-01_massHNL12p0_Vall2p314e-03-2016/nano_*_Friend.root")
#rdf = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/12Feb21/2016/HNL_majorana_all_ctau1p0e-01_massHNL12p0_Vall2p314e-03-2016/nano_*_Friend.root")

rdf_mu_coupling = rdf.Filter("LHEWeights_coupling_12>0")
rdf_e_coupling = rdf.Filter("LHEWeights_coupling_2>0")
rdf_tau_coupling = rdf.Filter("LHEWeights_coupling_67>0")

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_pt_ell_1_mu_coupling = rdf_mu_coupling.Histo1D(("hist_pt_ell_1_mu_coupling", "hist_pt_ell_1_mu_coupling", 25, 0., 100.), "GenLeptonsFromV_pt")
hist_pt_ell_1_e_coupling = rdf_e_coupling.Histo1D(("hist_pt_ell_1_e_coupling", "hist_pt_ell_1_e_coupling", 25, 0., 100.), "GenLeptonsFromV_pt")
hist_pt_ell_1_tau_coupling = rdf_tau_coupling.Histo1D(("hist_pt_ell_1_tau_coupling", "hist_pt_ell_1_tau_coupling", 24, 0., 100.), "GenLeptonsFromV_pt")

hist_pt_ell_1_e_coupling.SetLineStyle(2)
hist_pt_ell_1_tau_coupling.SetLineStyle(3)

hist_pt_ell_1_mu_coupling.SetLineWidth(2)
hist_pt_ell_1_e_coupling.SetLineWidth(2)
hist_pt_ell_1_tau_coupling.SetLineWidth(2)

hist_pt_ell_1_mu_coupling.Scale(1./hist_pt_ell_1_mu_coupling.Integral())
hist_pt_ell_1_e_coupling.Scale(1./hist_pt_ell_1_e_coupling.Integral())
hist_pt_ell_1_tau_coupling.Scale(1./hist_pt_ell_1_tau_coupling.Integral())

hist_pt_ell_1_mu_coupling.SetLineColor(color1)
hist_pt_ell_1_e_coupling.SetLineColor(color2)
hist_pt_ell_1_tau_coupling.SetLineColor(color3)

hist_pt_ell_1_mu_coupling.GetXaxis().SetTitle("p_{T} (l_{1}) [GeV]")
hist_pt_ell_1_mu_coupling.GetYaxis().SetTitle("Normalised entries")

hist_pt_ell_1_mu_coupling.GetYaxis().SetRangeUser(0, 1.2*max(hist_pt_ell_1_mu_coupling.GetMaximum(), hist_pt_ell_1_e_coupling.GetMaximum(), hist_pt_ell_1_tau_coupling.GetMaximum()))


hist_pt_ell_1_mu_coupling.Draw("HIST")
hist_pt_ell_1_e_coupling.Draw("HIST SAME")
hist_pt_ell_1_tau_coupling.Draw("HIST SAME")

leg = style.makeLegend(0.60, 0.6, 0.80, 0.78)
leg.AddEntry("hist_pt_ell_1_e_coupling", "e-only coupling", "l")
leg.AddEntry("hist_pt_ell_1_mu_coupling", "#mu-only coupling", "l")
leg.AddEntry("hist_pt_ell_1_tau_coupling", "#tau-only coupling", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.6, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeText(0.3, 0.78, 0.6, 0.8, "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm")

cv.SaveAs("kinematics/l1_pt.pdf")
#==========#
#==========#
#==========#

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_eta_ell_1_mu_coupling = rdf_mu_coupling.Histo1D(("hist_eta_ell_1_mu_coupling", "hist_eta_ell_1_mu_coupling", 25, -5, 5), "GenLeptonsFromV_eta")
hist_eta_ell_1_e_coupling = rdf_e_coupling.Histo1D(("hist_eta_ell_1_e_coupling", "hist_eta_ell_1_e_coupling", 25, -5, 5), "GenLeptonsFromV_eta")
hist_eta_ell_1_tau_coupling = rdf_tau_coupling.Histo1D(("hist_eta_ell_1_tau_coupling", "hist_eta_ell_1_tau_coupling", 24, -5, 5), "GenLeptonsFromV_eta")

hist_eta_ell_1_e_coupling.SetLineStyle(2)
hist_eta_ell_1_tau_coupling.SetLineStyle(3)

hist_eta_ell_1_mu_coupling.SetLineWidth(2)
hist_eta_ell_1_e_coupling.SetLineWidth(2)
hist_eta_ell_1_tau_coupling.SetLineWidth(2)

hist_eta_ell_1_mu_coupling.Scale(1./hist_eta_ell_1_mu_coupling.Integral())
hist_eta_ell_1_e_coupling.Scale(1./hist_eta_ell_1_e_coupling.Integral())
hist_eta_ell_1_tau_coupling.Scale(1./hist_eta_ell_1_tau_coupling.Integral())

hist_eta_ell_1_mu_coupling.SetLineColor(color1)
hist_eta_ell_1_e_coupling.SetLineColor(color2)
hist_eta_ell_1_tau_coupling.SetLineColor(color3)

hist_eta_ell_1_mu_coupling.GetXaxis().SetTitle("#eta (l_{1}) [GeV]")
hist_eta_ell_1_mu_coupling.GetYaxis().SetTitle("Normalised entries")

hist_eta_ell_1_mu_coupling.GetYaxis().SetRangeUser(0, 1.6*max(hist_eta_ell_1_mu_coupling.GetMaximum(), hist_eta_ell_1_e_coupling.GetMaximum(), hist_eta_ell_1_tau_coupling.GetMaximum()))


hist_eta_ell_1_mu_coupling.Draw("HIST")
hist_eta_ell_1_e_coupling.Draw("HIST SAME")
hist_eta_ell_1_tau_coupling.Draw("HIST SAME")

leg = style.makeLegend(0.60, 0.6, 0.80, 0.78)
leg.AddEntry("hist_eta_ell_1_e_coupling", "e-only coupling", "l")
leg.AddEntry("hist_eta_ell_1_mu_coupling", "#mu-only coupling", "l")
leg.AddEntry("hist_eta_ell_1_tau_coupling", "#tau-only coupling", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.6, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeText(0.3, 0.78, 0.6, 0.8, "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm")

cv.SaveAs("kinematics/l1_eta.pdf")

#==========#
#==========#
#==========#

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_pt_ell_2_mu_coupling = rdf_mu_coupling.Histo1D(("hist_pt_ell_2_mu_coupling", "hist_pt_ell_2_mu_coupling", 25, 0., 50.), "lepton_pt")
hist_pt_ell_2_e_coupling = rdf_e_coupling.Histo1D(("hist_pt_ell_2_e_coupling", "hist_pt_ell_2_e_coupling", 25, 0., 50.), "lepton_pt")
hist_pt_ell_2_tau_coupling = rdf_tau_coupling.Histo1D(("hist_pt_ell_2_tau_coupling", "hist_pt_ell_2_tau_coupling", 24, 0., 50.), "lepton_pt")

hist_pt_ell_2_e_coupling.SetLineStyle(2)
hist_pt_ell_2_tau_coupling.SetLineStyle(3)

hist_pt_ell_2_mu_coupling.SetLineWidth(2)
hist_pt_ell_2_e_coupling.SetLineWidth(2)
hist_pt_ell_2_tau_coupling.SetLineWidth(2)

hist_pt_ell_2_mu_coupling.Scale(1./hist_pt_ell_2_mu_coupling.Integral())
hist_pt_ell_2_e_coupling.Scale(1./hist_pt_ell_2_e_coupling.Integral())
hist_pt_ell_2_tau_coupling.Scale(1./hist_pt_ell_2_tau_coupling.Integral())

hist_pt_ell_2_mu_coupling.SetLineColor(color1)
hist_pt_ell_2_e_coupling.SetLineColor(color2)
hist_pt_ell_2_tau_coupling.SetLineColor(color3)

hist_pt_ell_2_mu_coupling.GetXaxis().SetTitle("p_{T} (l_{2}) [GeV]")
hist_pt_ell_2_mu_coupling.GetYaxis().SetTitle("Normalised entries")

hist_pt_ell_2_mu_coupling.GetYaxis().SetRangeUser(0, 1.2*max(hist_pt_ell_2_mu_coupling.GetMaximum(), hist_pt_ell_2_e_coupling.GetMaximum(), hist_pt_ell_2_tau_coupling.GetMaximum()))


hist_pt_ell_2_mu_coupling.Draw("HIST")
hist_pt_ell_2_e_coupling.Draw("HIST SAME")
hist_pt_ell_2_tau_coupling.Draw("HIST SAME")

leg = style.makeLegend(0.60, 0.6, 0.80, 0.78)
leg.AddEntry("hist_pt_ell_2_e_coupling", "e-only coupling", "l")
leg.AddEntry("hist_pt_ell_2_mu_coupling", "#mu-only coupling", "l")
leg.AddEntry("hist_pt_ell_2_tau_coupling", "#tau-only coupling", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.6, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeText(0.3, 0.78, 0.6, 0.8, "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm")

cv.SaveAs("kinematics/l2_pt.pdf")
#==========#
#==========#
#==========#

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_eta_ell_2_mu_coupling = rdf_mu_coupling.Histo1D(("hist_eta_ell_2_mu_coupling", "hist_eta_ell_2_mu_coupling", 25, -5, 5), "lepton_eta")
hist_eta_ell_2_e_coupling = rdf_e_coupling.Histo1D(("hist_eta_ell_2_e_coupling", "hist_eta_ell_2_e_coupling", 25, -5, 5), "lepton_eta")
hist_eta_ell_2_tau_coupling = rdf_tau_coupling.Histo1D(("hist_eta_ell_2_tau_coupling", "hist_eta_ell_2_tau_coupling", 24, -5, 5), "lepton_eta")

hist_eta_ell_2_e_coupling.SetLineStyle(2)
hist_eta_ell_2_tau_coupling.SetLineStyle(3)

hist_eta_ell_2_mu_coupling.SetLineWidth(2)
hist_eta_ell_2_e_coupling.SetLineWidth(2)
hist_eta_ell_2_tau_coupling.SetLineWidth(2)

hist_eta_ell_2_mu_coupling.Scale(1./hist_eta_ell_2_mu_coupling.Integral())
hist_eta_ell_2_e_coupling.Scale(1./hist_eta_ell_2_e_coupling.Integral())
hist_eta_ell_2_tau_coupling.Scale(1./hist_eta_ell_2_tau_coupling.Integral())

hist_eta_ell_2_mu_coupling.SetLineColor(color1)
hist_eta_ell_2_e_coupling.SetLineColor(color2)
hist_eta_ell_2_tau_coupling.SetLineColor(color3)

hist_eta_ell_2_mu_coupling.GetXaxis().SetTitle("#eta (l_{2}) [GeV]")
hist_eta_ell_2_mu_coupling.GetYaxis().SetTitle("Normalised entries")

hist_eta_ell_2_mu_coupling.GetYaxis().SetRangeUser(0, 1.6*max(hist_eta_ell_2_mu_coupling.GetMaximum(), hist_eta_ell_2_e_coupling.GetMaximum(), hist_eta_ell_2_tau_coupling.GetMaximum()))


hist_eta_ell_2_mu_coupling.Draw("HIST")
hist_eta_ell_2_e_coupling.Draw("HIST SAME")
hist_eta_ell_2_tau_coupling.Draw("HIST SAME")

leg = style.makeLegend(0.60, 0.6, 0.80, 0.78)
leg.AddEntry("hist_eta_ell_2_e_coupling", "e-only coupling", "l")
leg.AddEntry("hist_eta_ell_2_mu_coupling", "#mu-only coupling", "l")
leg.AddEntry("hist_eta_ell_2_tau_coupling", "#tau-only coupling", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.6, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
style.makeText(0.3, 0.78, 0.6, 0.8, "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm")

cv.SaveAs("kinematics/l2_eta.pdf")

#==========#
#==========#
#==========#
rdf_1 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/leptonEff/03Nov20/HNL_majorana_all_ctau1p0e02_massHNL1p0_Vall5p274e-02-2016/nano_*_Friend.root")
rdf_6 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/leptonEff/03Nov20/HNL_majorana_all_ctau1p0e01_massHNL6p0_Vall1p454e-03-2016/nano_*_Friend.root")
rdf_12 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/leptonEff/03Nov20/HNL_majorana_all_ctau1p0e-01_massHNL12p0_Vall2p314e-03-2016/nano_*_Friend.root")


#rdf_1 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/12Feb21/HNL_majorana_all_ctau1p0e02_massHNL1p0_Vall5p274e-02-2016/2016/nano_*_Friend.root")
#rdf_6 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/12Feb21/HNL_majorana_all_ctau1p0e01_massHNL6p0_Vall1p454e-03-2016/2016/nano_*_Friend.root")
#rdf_12 = ROOT.RDataFrame("Friends", "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/12Feb21/HNL_majorana_all_ctau1p0e-01_massHNL12p0_Vall2p314e-03-2016/2016/nano_*_Friend.root")

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
rdf_mu_coupling_filtered = rdf_mu_coupling.Filter("nGenLeptonsFromV==1").Filter("nlepton==1")
rdf_mu_coupling_filtered = rdf_mu_coupling_filtered.Define("l1_pt", "GenLeptonsFromV_pt[0]")
rdf_mu_coupling_filtered = rdf_mu_coupling_filtered.Define("l2_pt", "lepton_pt[0]")

#norm = rdf_mu_coupling_filtered.Filter("l1_pt>l2_pt").Count().GetValue()/rdf_mu_coupling_filtered.Count().GetValue()
#print(norm)

hist_pt_ell_1_pt_ell_2_mu_coupling = rdf_mu_coupling_filtered.Histo2D(("hist_pt_ell_1_pt_ell_2_mu_coupling", "hist_pt_ell_1_pt_ell_2_mu_coupling", 25, 0, 50, 25, 0, 50 ), "l1_pt", "l2_pt")

hist_pt_ell_1_pt_ell_2_mu_coupling.GetXaxis().SetTitle("p_{T} (l_{1}) [GeV]")
hist_pt_ell_1_pt_ell_2_mu_coupling.GetYaxis().SetTitle("p_{T} (l_{2}) [GeV]")

hist_pt_ell_1_pt_ell_2_mu_coupling.Draw("COLZ")
l = ROOT.TLine(0., 0., 50., 50.)
l.SetLineColor(ROOT.kWhite)
l.SetLineStyle(2)
l.SetLineWidth(3)
l.Draw("")

#style.makeText(0.3, 0.83, 0.6, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightar3owNl_{1}")6
cv.SaveAs("kinematics/l1_pt_l2_pt.pdf")

#==========#
#==========#
cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_pt_HNL_1 = rdf_1.Histo1D(("hist_pt_HNL_1", "hist_pt_HNL_1", 25, 0., 50.), "HNL_pt")
hist_pt_HNL_6 = rdf_6.Histo1D(("hist_pt_HNL_6", "hist_pt_HNL_6", 25, 0., 50.), "HNL_pt")
hist_pt_HNL_12= rdf_12.Histo1D(("hist_pt_HNL_12", "hist_pt_HNL_12", 25, 0., 50.), "HNL_pt")

hist_pt_HNL_6.SetLineStyle(2)
hist_pt_HNL_12.SetLineStyle(3)

hist_pt_HNL_1.SetLineWidth(2)
hist_pt_HNL_6.SetLineWidth(2)
hist_pt_HNL_12.SetLineWidth(2)

hist_pt_HNL_1.Scale(1./hist_pt_HNL_1.Integral())
hist_pt_HNL_6.Scale(1./hist_pt_HNL_6.Integral())
hist_pt_HNL_12.Scale(1./hist_pt_HNL_12.Integral())

hist_pt_HNL_1.SetLineColor(color1)
hist_pt_HNL_6.SetLineColor(color2)
hist_pt_HNL_12.SetLineColor(color3)

hist_pt_HNL_1.GetXaxis().SetTitle("p_{T} (N) [GeV]")
hist_pt_HNL_1.GetYaxis().SetTitle("Normalised entries")

hist_pt_HNL_1.GetYaxis().SetRangeUser(0, 1.6*max(hist_pt_HNL_1.GetMaximum(), hist_pt_HNL_6.GetMaximum(), hist_pt_HNL_12.GetMaximum()))

hist_pt_HNL_1.Draw("HIST")
hist_pt_HNL_6.Draw("HIST SAME")
hist_pt_HNL_12.Draw("HIST SAME")

leg = style.makeLegend(0.3, 0.65, 0.7, 0.82)
leg.AddEntry("hist_pt_HNL_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_pt_HNL_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_pt_HNL_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.6, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightar3owNl_{1}")
cv.SaveAs("kinematics/HNL_pt.pdf")



#==========#
#==========#
cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
hist_eta_HNL_1 = rdf_1.Histo1D(("hist_eta_HNL_1", "hist_eta_HNL_1", 25, -5, 5), "HNL_eta")
hist_eta_HNL_6 = rdf_6.Histo1D(("hist_eta_HNL_6", "hist_eta_HNL_6", 25, -5, 5), "HNL_eta")
hist_eta_HNL_12= rdf_12.Histo1D(("hist_eta_HNL_12", "hist_eta_HNL_12", 25, -5, 5), "HNL_eta")

hist_eta_HNL_6.SetLineStyle(2)
hist_eta_HNL_12.SetLineStyle(3)

hist_eta_HNL_1.SetLineWidth(2)
hist_eta_HNL_6.SetLineWidth(2)
hist_eta_HNL_12.SetLineWidth(2)

hist_eta_HNL_1.Scale(1./hist_eta_HNL_1.Integral())
hist_eta_HNL_6.Scale(1./hist_eta_HNL_6.Integral())
hist_eta_HNL_12.Scale(1./hist_eta_HNL_12.Integral())

hist_eta_HNL_1.SetLineColor(color1)
hist_eta_HNL_6.SetLineColor(color2)
hist_eta_HNL_12.SetLineColor(color3)

hist_eta_HNL_1.GetXaxis().SetTitle("#eta (N) [GeV]")
hist_eta_HNL_1.GetYaxis().SetTitle("Normalised entries")

hist_eta_HNL_1.GetYaxis().SetRangeUser(0, 1.6*max(hist_eta_HNL_1.GetMaximum(), hist_eta_HNL_6.GetMaximum(), hist_eta_HNL_12.GetMaximum()))

hist_eta_HNL_1.Draw("HIST")
hist_eta_HNL_6.Draw("HIST SAME")
hist_eta_HNL_12.Draw("HIST SAME")

leg = style.makeLegend(0.3, 0.65, 0.7, 0.82)
leg.AddEntry("hist_eta_HNL_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_eta_HNL_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_eta_HNL_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")

style.makeText(0.3, 0.83, 0.5, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
cv.SaveAs("kinematics/HNL_eta.pdf")

#==========#
#==========#
cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
cv.SetLogx()
bins = np.geomspace(1, 100, num=26)
hist_boost_HNL_1 = rdf_1.Histo1D(("hist_boost_HNL_1", "hist_boost_HNL_1", 25, bins), "HNL_boost")
hist_boost_HNL_6 = rdf_6.Histo1D(("hist_boost_HNL_6", "hist_boost_HNL_6", 25, bins), "HNL_boost")
hist_boost_HNL_12= rdf_12.Histo1D(("hist_boost_HNL_12", "hist_boost_HNL_12", 25, bins), "HNL_boost")

hist_boost_HNL_6.SetLineStyle(2)
hist_boost_HNL_12.SetLineStyle(3)

hist_boost_HNL_1.SetLineWidth(2)
hist_boost_HNL_6.SetLineWidth(2)
hist_boost_HNL_12.SetLineWidth(2)

hist_boost_HNL_1.Scale(1./hist_boost_HNL_1.Integral())
hist_boost_HNL_6.Scale(1./hist_boost_HNL_6.Integral())
hist_boost_HNL_12.Scale(1./hist_boost_HNL_12.Integral())

hist_boost_HNL_1.SetLineColor(color1)
hist_boost_HNL_6.SetLineColor(color2)
hist_boost_HNL_12.SetLineColor(color3)

hist_boost_HNL_1.GetXaxis().SetTitle("#beta#gamma (N)")
hist_boost_HNL_1.GetYaxis().SetTitle("Normalised entries")

hist_boost_HNL_1.GetYaxis().SetRangeUser(0, 1.6*max(hist_boost_HNL_1.GetMaximum(), hist_boost_HNL_6.GetMaximum(), hist_boost_HNL_12.GetMaximum()))

hist_boost_HNL_1.Draw("HIST")
hist_boost_HNL_6.Draw("HIST SAME")
hist_boost_HNL_12.Draw("HIST SAME")

leg = style.makeLegend(0.3, 0.65, 0.7, 0.82)
leg.AddEntry("hist_boost_HNL_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_boost_HNL_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_boost_HNL_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.5, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")
cv.SaveAs("kinematics/HNL_boost.pdf")

#==========#
#==========#
cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.15)

cv.Draw("")
cv.SetLogx()

bins = np.geomspace(0.0001, 1000000, num=51)
hist_dxy_HNL_1 = rdf_1.Histo1D(("hist_dxy_HNL_1", "hist_dxy_HNL_1", 50, bins), "HNL_dxy")
hist_dxy_HNL_6 = rdf_6.Histo1D(("hist_dxy_HNL_6", "hist_dxy_HNL_6", 50, bins), "HNL_dxy")
hist_dxy_HNL_12= rdf_12.Histo1D(("hist_dxy_HNL_12", "hist_dxy_HNL_12", 50, bins), "HNL_dxy")

hist_dxy_HNL_6.SetLineStyle(2)
hist_dxy_HNL_12.SetLineStyle(3)

hist_dxy_HNL_1.SetLineWidth(2)
hist_dxy_HNL_6.SetLineWidth(2)
hist_dxy_HNL_12.SetLineWidth(2)

hist_dxy_HNL_1.Scale(1./hist_dxy_HNL_1.Integral())
hist_dxy_HNL_6.Scale(1./hist_dxy_HNL_6.Integral())
hist_dxy_HNL_12.Scale(1./hist_dxy_HNL_12.Integral())

hist_dxy_HNL_1.SetLineColor(color1)
hist_dxy_HNL_6.SetLineColor(color2)
hist_dxy_HNL_12.SetLineColor(color3)

hist_dxy_HNL_1.GetXaxis().SetTitle("L_{xy} (N) [cm]")
hist_dxy_HNL_1.GetYaxis().SetTitle("Normalised entries")

hist_dxy_HNL_1.GetYaxis().SetRangeUser(0, 1.6*max(hist_dxy_HNL_1.GetMaximum(), hist_dxy_HNL_6.GetMaximum(), hist_dxy_HNL_12.GetMaximum()))

hist_dxy_HNL_1.Draw("HIST")
hist_dxy_HNL_6.Draw("HIST SAME")
hist_dxy_HNL_12.Draw("HIST SAME")

leg = style.makeLegend(0.3, 0.65, 0.7, 0.82)
leg.AddEntry("hist_dxy_HNL_1", "m_{N} = 1 GeV, c#tau_{0} = 10 cm", "l")
leg.AddEntry("hist_dxy_HNL_6", "m_{N} = 6 GeV, c#tau_{0} = 1 cm", "l")
leg.AddEntry("hist_dxy_HNL_12", "m_{N} = 12 GeV, c#tau_{0} = 0.1 mm", "l")
leg.Draw("")
style.makeText(0.3, 0.83, 0.5, 0.83, "#sqrt{s}=13 TeV, pp#rightarrowW#rightarrowNl_{1}")

cv.SaveAs("kinematics/HNL_dxy.pdf")













