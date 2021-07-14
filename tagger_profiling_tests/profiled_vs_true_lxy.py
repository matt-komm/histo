import ROOT
import os 
import numpy as np
from histo import style

file_list = ROOT.std.vector('string')()
path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/28May21/2016"
for directory in os.listdir(path):
    if "pt20" in directory:
        for f in os.listdir(os.path.join(path, directory)):
            file_list.push_back(os.path.join(path, directory, f))

rdf = ROOT.RDataFrame("Friends", file_list)
rdf = rdf.Define("hnlJet_nominal_truelogLxy", "log10(hnlJet_nominal_trueLxy)")
rdf = rdf.Filter("nselectedJets_nominal>0")
bin_ranges = np.linspace(-2., 2, num=20)
bins = ("", "", 19, bin_ranges, 19, bin_ranges)
for jet_class in ["Q", "QMU", "QE"]:
    hist = rdf.Histo2D(bins, "hnlJet_nominal_trueLxy", f"hnlJet_nominal_llpdnnx_ratio_LLP_{jet_class}_param", f"hnlJet_nominal_isTrue{jet_class}")
    cv = style.makeCanvas("")
    cv.Draw("")
    cv.SetBottomMargin(0.15)
    cv.SetLeftMargin(0.10)
    cv.SetRightMargin(0.15)

    hist.GetXaxis().SetTitle("Log_{10} (True L_{xy}/1cm)")
    hist.GetYaxis().SetTitle("Log_{10} (Profiled #hat{L}_{xy}/1cm)")
    hist.GetYaxis().SetTitleOffset(0.6*hist.GetYaxis().GetTitleOffset())

    hist.Scale(1./hist.Integral())
    hist.Draw("COLZ")
    l = ROOT.TLine(-2, -2, 2, 2)
    l.SetLineStyle(3)
    l.SetLineWidth(2)
    l.Draw("SAME")
    style.makeText(0.2, 0.7, 0.38, 0.8, f"LLP_{jet_class}")
    cv.SaveAs(f"true_vs_profile_lxy_{jet_class}.pdf")
    cv.SaveAs(f"true_vs_profile_lxy_{jet_class}.png")
