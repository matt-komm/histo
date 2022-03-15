import ROOT
import numpy as np
from histo import style
import os
import matplotlib.pyplot as plt

color1 = ROOT.kOrange+7
color2 = ROOT.kSpring+3

file_list = ROOT.std.vector('string')()
path = "/vols/cms/vc1117/LLP/nanoAOD_friends/leptonEff/11Sep21"
for directory in os.listdir(path):
    for f in os.listdir(os.path.join(path, directory)):
        file_list.push_back(os.path.join(path, directory, f))

rdf = ROOT.RDataFrame("Friends", file_list)

rdf = rdf.Define("muon_jetmatch", "(muon_pt>3.)&&abs(muon_eta)<2.4")
rdf = rdf.Define("muon_match", "muon_jetmatch*muon_is_matched&&(muon_reco_pt>3.)&&abs(muon_reco_eta)<2.4")
rdf = rdf.Define("muon_id", "muon_match*muon_looseId")

rdf = rdf.Define("electron_jetmatch", "(electron_pt>5.)&&abs(electron_eta)<2.4")
rdf = rdf.Define("electron_match", "electron_jetmatch*electron_is_matched&&(electron_reco_pt>5.)&&abs(electron_reco_eta)<2.4")
rdf = rdf.Define("electron_id", "electron_match*electron_customID")

rdf = rdf.Define("jet_jetmatch", "(jet_pt>15.) && abs(jet_eta)<2.4")
rdf = rdf.Define("jet_match", "jet_jetmatch*jet_is_matched&&(jet_reco_pt>15.)&&abs(jet_reco_eta)<2.4")
rdf = rdf.Define("jet_id", "jet_match*jet_tightId")

bin_ranges = np.logspace(-1., 3, num=30)
bins = ("", "", 29, bin_ranges)
'''

mu_numerator = rdf.Histo1D(bins, "muon_dxy", "muon_match")
mu_denominator = rdf.Histo1D(bins, "muon_dxy", "muon_jetmatch")
mu_id_numerator = rdf.Histo1D(bins, "muon_dxy", "muon_id")
mu_id_denominator = rdf.Histo1D(bins, "muon_dxy", "muon_match")

el_numerator = rdf.Histo1D(bins, "electron_dxy", "electron_match")
el_denominator = rdf.Histo1D(bins, "electron_dxy", "electron_jetmatch")
el_id_numerator = rdf.Histo1D(bins, "electron_dxy", "electron_id")
el_id_denominator = rdf.Histo1D(bins, "electron_dxy", "electron_match")

jet_numerator = rdf.Histo1D(bins, "jet_dxy", "jet_match")
jet_denominator = rdf.Histo1D(bins, "jet_dxy", "jet_jetmatch")
jet_id_numerator = rdf.Histo1D(bins, "jet_dxy", "jet_id")
jet_id_denominator = rdf.Histo1D(bins, "jet_dxy", "jet_match")


mu_eff = ROOT.TEfficiency(mu_numerator.Clone(), mu_denominator.Clone())
el_eff = ROOT.TEfficiency(el_numerator.Clone(), el_denominator.Clone())
jet_eff = ROOT.TEfficiency(jet_numerator.Clone(), jet_denominator.Clone())
mu_id_eff = ROOT.TEfficiency(mu_id_numerator.Clone(), mu_id_denominator.Clone())
el_id_eff = ROOT.TEfficiency(el_id_numerator.Clone(), el_id_denominator.Clone())
jet_id_eff = ROOT.TEfficiency(jet_id_numerator.Clone(), jet_id_denominator.Clone())

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
jet_eff.Draw("A")
jet_id_eff.Draw("SAME")
jet_eff.SetLineColor(color1)
jet_eff.SetMarkerColor(color1)
jet_id_eff.SetLineColor(color2)
jet_id_eff.SetMarkerColor(color2)
jet_eff.SetLineWidth(3)
jet_id_eff.SetLineStyle(2)
jet_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.2, 0.25, 0.35, 0.35)
leg.SetTextSize(24)
leg.AddEntry(jet_eff, "jet reco/gen", "pl")
leg.AddEntry(jet_id_eff, "jet id/reco", "pl")
leg.Draw("SAME")

cv.Update()
cv.SetLogx()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.2, 0.45, 0.2, 0.45, "p_{T} > 15 GeV, |#eta| < 2.4 (gen&reco)", size=24)
style.makeText(0.2, 0.5, 0.2, 0.5, "#sqrt{s}=13 TeV, N#rightarrow#mujj", size=24)

jet_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
jet_eff.GetPaintedGraph().GetXaxis().SetRangeUser(0., 300.)
jet_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)

jet_eff.SetTitle(";L^{gen}_{xy} [cm];Efficiency")

cv.SaveAs("jet_dxy.pdf")


cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
mu_eff.Draw("A")
mu_id_eff.Draw("SAME")
mu_eff.SetLineColor(color1)
mu_eff.SetMarkerColor(color1)
mu_id_eff.SetLineColor(color2)
mu_id_eff.SetMarkerColor(color2)
mu_eff.SetLineWidth(3)
mu_id_eff.SetLineStyle(2)
mu_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.2, 0.25, 0.35, 0.35)
leg.SetTextSize(24)
leg.AddEntry(mu_eff, "#mu reco/gen", "pl")
leg.AddEntry(mu_id_eff, "#mu id/reco", "pl")
leg.Draw("SAME")

cv.Update()
cv.SetLogx()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.2, 0.45, 0.2, 0.45, "p_{T} > 3 GeV, |#eta| < 2.4 (gen&reco)", size=24)
style.makeText(0.2, 0.5, 0.2, 0.5, "#sqrt{s}=13 TeV, N#rightarrow#mujj", size=24)

mu_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
mu_eff.GetPaintedGraph().GetXaxis().SetRangeUser(0., 100.)
mu_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)

mu_eff.SetTitle(";L^{gen}_{xy} [cm];Efficiency")

cv.SaveAs("mu_dxy.pdf")

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
el_eff.Draw("A")
el_id_eff.Draw("SAME")
el_eff.SetLineColor(color1)
el_eff.SetMarkerColor(color1)
el_id_eff.SetLineColor(color2)
el_id_eff.SetMarkerColor(color2)
el_eff.SetLineWidth(3)
el_id_eff.SetLineStyle(2)
el_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.17, 0.25, 0.32, 0.35)
leg.SetTextSize(24)
leg.AddEntry(el_eff, "e reco/gen", "pl")
leg.AddEntry(el_id_eff, "e id/reco", "pl")
leg.Draw("SAME")

cv.Update()
cv.SetLogx()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.17, 0.45, 0.17, 0.45, "p_{T} > 5 GeV, |#eta| < 2.4 (gen&reco)", size=24)
style.makeText(0.17, 0.5, 0.17, 0.5, "#sqrt{s}=13 TeV, N#rightarrowejj", size=24)

el_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
el_eff.GetPaintedGraph().GetXaxis().SetRangeUser(0., 100.)
el_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)

el_eff.SetTitle(";L^{gen}_{xy} [cm];Efficiency")

cv.SaveAs("el_dxy.pdf")
'''

#### PT #####

rdf = ROOT.RDataFrame("Friends", file_list)

rdf = rdf.Define("muon_jetmatch", "(muon_pt>3.)&&abs(muon_eta)<2.4&&(muon_dxy<1)")
rdf = rdf.Define("muon_match", "muon_jetmatch*muon_is_matched&&abs(muon_reco_eta)<2.4")
rdf = rdf.Define("muon_id", "muon_match*muon_looseId")

rdf = rdf.Define("electron_jetmatch", "(electron_pt>3.)&&abs(electron_eta)<2.4&&(electron_dxy<1)")
rdf = rdf.Define("electron_match", "electron_jetmatch*electron_is_matched&&abs(electron_reco_eta)<2.4")
rdf = rdf.Define("electron_id", "electron_match*electron_customID")

rdf = rdf.Define("jet_jetmatch", "(jet_pt>15.)&&abs(jet_eta)<2.4")
rdf = rdf.Define("jet_match", "jet_jetmatch*jet_is_matched&&abs(jet_reco_eta)<2.4")
rdf = rdf.Define("jet_id", "jet_match*jet_tightId")

rdf = rdf.Define("cpf_weighted", "jet_match*jet_cpf_sum")
rdf = rdf.Define("npf_weighted", "jet_match*jet_npf_sum")

bin_ranges_mu = np.linspace(3., 13., num=20)
bins_mu = ("", "", 19, bin_ranges_mu)

bin_ranges_e = np.linspace(5., 15., num=20)
bins_e = ("", "", 19, bin_ranges_e)

bin_ranges_jet = np.linspace(15., 60., num=20)
bins_jet = ("", "", 19, bin_ranges_jet)

bins_energy = ("", "", 19, np.linspace(0., 1., num=20))

mu_numerator = rdf.Histo1D(bins_mu, "muon_pt", "muon_match")
mu_denominator = rdf.Histo1D(bins_mu, "muon_pt", "muon_jetmatch")
mu_id_numerator = rdf.Histo1D(bins_mu, "muon_pt", "muon_id")
mu_id_denominator = rdf.Histo1D(bins_mu, "muon_pt", "muon_match")

el_numerator = rdf.Histo1D(bins_e, "electron_pt", "electron_match")
el_denominator = rdf.Histo1D(bins_e, "electron_pt", "electron_jetmatch")
el_id_numerator = rdf.Histo1D(bins_e, "electron_pt", "electron_id")
el_id_denominator = rdf.Histo1D(bins_e, "electron_pt", "electron_match")

jet_numerator = rdf.Histo1D(bins_jet, "jet_pt", "jet_match")
jet_denominator = rdf.Histo1D(bins_jet, "jet_pt", "jet_jetmatch")
jet_id_numerator = rdf.Histo1D(bins_jet, "jet_pt", "jet_id")
jet_id_denominator = rdf.Histo1D(bins_jet, "jet_pt", "jet_match")

jet_cpf = rdf.Histo1D(bins, "jet_dxy", "cpf_weighted")
jet_npf = rdf.Histo1D(bins, "jet_dxy", "npf_weighted")
jet_denominator = rdf.Histo1D(bins, "jet_dxy", "jet_match")

jet_cpf.Divide(jet_denominator.Clone())
jet_npf.Divide(jet_denominator.Clone())

jet_cpf = jet_cpf.Clone("cpf")
jet_npf = jet_npf.Clone("npf")


cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
jet_cpf.Draw("")
jet_npf.Draw("SAME")
jet_cpf.SetLineColor(color1)
jet_cpf.SetMarkerColor(color1)
jet_npf.SetLineColor(color2)
jet_npf.SetMarkerColor(color2)

leg = style.makeLegend(0.2, 0.2, 0.35, 0.3)
leg.SetTextSize(24)
leg.AddEntry(jet_cpf, "charged particle", "pl")
leg.AddEntry(jet_npf, "neutral particle", "pl")
leg.Draw("SAME")

cv.Update()
cv.SetLogx()
#style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.2, 0.75, 0.2, 0.75, "p_{T} > 15 GeV, |#eta| < 2.4 (gen&reco)", size=24)
style.makeText(0.2, 0.8, 0.2, 0.8, "#sqrt{s}=13 TeV, N#rightarrow#mujj", size=24)

jet_cpf.GetYaxis().SetRangeUser(0, 1.)
jet_cpf.GetXaxis().SetRangeUser(0., 300.)
jet_cpf.GetXaxis().SetNdivisions(510)

jet_cpf.SetTitle(";L^{gen}_{xy} [cm];Average energy fraction")

cv.SaveAs("jet_energy_fraction.pdf")

'''
mu_eff = ROOT.TEfficiency(mu_numerator.Clone(), mu_denominator.Clone())
el_eff = ROOT.TEfficiency(el_numerator.Clone(), el_denominator.Clone())
jet_eff = ROOT.TEfficiency(jet_numerator.Clone(), jet_denominator.Clone())

mu_id_eff = ROOT.TEfficiency(mu_id_numerator.Clone(), mu_id_denominator.Clone())
el_id_eff = ROOT.TEfficiency(el_id_numerator.Clone(), el_id_denominator.Clone())
jet_id_eff = ROOT.TEfficiency(jet_id_numerator.Clone(), jet_id_denominator.Clone())

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
jet_eff.Draw("A")
jet_id_eff.Draw("SAME")
jet_eff.SetLineColor(color1)
jet_eff.SetMarkerColor(color1)
jet_id_eff.SetLineColor(color2)
jet_id_eff.SetMarkerColor(color2)
jet_eff.SetLineWidth(3)
jet_id_eff.SetLineStyle(2)
jet_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.4, 0.25, 0.55, 0.35)
leg.SetTextSize(24)
leg.AddEntry(jet_eff, "jet reco/gen", "pl")
leg.AddEntry(jet_id_eff, "jet id/reco", "pl")
leg.Draw("SAME")

cv.Update()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.4, 0.45, 0.4, 0.45, "p_{T} > 15 GeV, |#eta| < 2.4, L_{xy} < 1 cm (gen)", size=24)
style.makeText(0.4, 0.5, 0.4, 0.5, "#sqrt{s}=13 TeV, N#rightarrow#mujj", size=24)

jet_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
jet_eff.GetPaintedGraph().GetXaxis().SetRangeUser(15., 60.)
jet_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)
jet_eff.SetTitle(";p^{gen}_{T} [GeV];Efficiency")

cv.SaveAs("jet_pt.pdf")


cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
mu_eff.Draw("A")
mu_id_eff.Draw("SAME")
mu_eff.SetLineColor(color1)
mu_eff.SetMarkerColor(color1)
mu_id_eff.SetLineColor(color2)
mu_id_eff.SetMarkerColor(color2)
mu_eff.SetLineWidth(3)
mu_id_eff.SetLineStyle(2)
mu_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.4, 0.25, 0.55, 0.35)
leg.SetTextSize(24)
leg.AddEntry(mu_eff, "#mu reco/gen", "pl")
leg.AddEntry(mu_id_eff, "#mu id/reco", "pl")
leg.Draw("SAME")

cv.Update()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.4, 0.45, 0.4, 0.45, "p_{T} > 3 GeV, |#eta| < 2.4, L_{xy} < 1 cm (gen)", size=24)
style.makeText(0.4, 0.5, 0.4, 0.5, "#sqrt{s}=13 TeV, N#rightarrow#mujj", size=24)

mu_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
mu_eff.GetPaintedGraph().GetXaxis().SetRangeUser(3., 13.)
mu_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)
mu_eff.SetTitle(";p^{gen}_{T} [GeV];Efficiency")

cv.SaveAs("mu_pt.pdf")


cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
mu_eff.Draw("A")
mu_id_eff.Draw("SAME")
mu_eff.SetLineColor(color1)
mu_eff.SetMarkerColor(color1)
mu_id_eff.SetLineColor(color2)
mu_id_eff.SetMarkerColor(color2)
mu_eff.SetLineWidth(3)
mu_id_eff.SetLineStyle(2)
mu_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.4, 0.25, 0.55, 0.35)
leg.SetTextSize(24)
leg.AddEntry(mu_eff, "#mu reco/gen", "pl")
leg.AddEntry(mu_id_eff, "#mu id/reco", "pl")
leg.Draw("SAME")

cv.Update()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.4, 0.45, 0.4, 0.45, "p_{T} > 3 GeV, |#eta| < 2.4, L_{xy} < 1 cm (gen)", size=24)
style.makeText(0.4, 0.5, 0.4, 0.5, "#sqrt{s}=13 TeV, N#rightarrow#mujj", size=24)

mu_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
mu_eff.GetPaintedGraph().GetXaxis().SetRangeUser(3., 13.)
mu_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)
mu_eff.SetTitle(";p^{gen}_{T} [GeV];Efficiency")

cv.SaveAs("mu_pt.pdf")

cv = style.makeCanvas()
cv.SetBottomMargin(0.15)
cv.SetLeftMargin(0.13)
cv.Draw("")
el_eff.Draw("A")
el_id_eff.Draw("SAME")
el_eff.SetLineColor(color1)
el_eff.SetMarkerColor(color1)
el_id_eff.SetLineColor(color2)
el_id_eff.SetMarkerColor(color2)
el_eff.SetLineWidth(3)
el_id_eff.SetLineStyle(2)
el_id_eff.SetLineWidth(3)

leg = style.makeLegend(0.42, 0.25, 0.77, 0.35)
leg.SetTextSize(24)
leg.AddEntry(el_eff, "e reco/gen", "pl")
leg.AddEntry(el_id_eff, "e id/reco", "pl")
leg.Draw("SAME")

cv.Update()
style.makeCMSText(0.15, 0.86, additionalText="Simulation Internal")
style.makeText(0.42, 0.45, 0.42, 0.45, "p_{T} > 3 GeV, |#eta| < 2.4, L_{xy} < 1 cm (gen)", size=24)
style.makeText(0.42, 0.5, 0.42, 0.5, "#sqrt{s}=13 TeV, N#rightarrowejj", size=24)

el_eff.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.2)
el_eff.GetPaintedGraph().GetXaxis().SetRangeUser(5., 15.)
el_eff.GetPaintedGraph().GetXaxis().SetNdivisions(510)

el_eff.SetTitle(";p^{gen}_{T} [GeV];Efficiency")

cv.SaveAs("el_pt.pdf")



'''