import ROOT
import os
from histo import style

color1 = ROOT.kOrange+7
color2 = ROOT.kSpring+3

benchmark_model = "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03-2016"
ntuple_path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/29Jun21/2016/" 

file_list = ROOT.std.vector('string')()

for f in os.listdir(os.path.join(ntuple_path, benchmark_model)):
    file_list.push_back(os.path.join(ntuple_path, benchmark_model, f))


uncertainty_list = {"pileup": "puweight_nominal",
                    "muon": "IsoMuTrigger_weight_trigger_nominal*tightMuons_weight_iso_nominal*tightMuons_weight_id_nominal",
                    "electron": "tightElectrons_weight_id_nominal*tightElectrons_weight_reco_nominal*looseElectrons_weight_reco_nominal",
                    "jesTotal": "shape",
                    "unclEn": "shape",
                    "jer": "shape",
                    "track_unc": "hnlJet_track_weight_adapted_nominal",
                    "pdf": "pdf_nominal",
                    "scale": "scale_shapeonly_nominal"
                    }

uncertainty_text = {"pileup": "pileup",
                    "muon": "muon combined",
                    "electron": "electron combined",
                    "jesTotal": "jet energy scale",
                    "unclEn": "unclustered energy",
                    "jer": "jet energy resolution",
                    "track_unc": "displaced track",
                    "pdf": "PDF",
                    "scale": "scale"
                    }

for final_state in ["mu", "e"]:
    rdf = ROOT.RDataFrame("Friends", file_list)
    rdf = rdf.Filter("nominal_dR_l2j<0.4")
    if final_state == "mu":
        rdf = rdf.Filter("subleadingLeptons_isMuon[0]")
        rdf = rdf.Define("discriminant", "hnlJet_nominal_llpdnnx_ratio_LLP_QMU")
    elif final_state == "e":
        rdf = rdf.Filter("subleadingLeptons_isElectron[0]")
        rdf = rdf.Define("discriminant", "hnlJet_nominal_llpdnnx_ratio_LLP_QE")
    for key, uncertainty in uncertainty_list.items():
        cv = style.makeCanvas()
        leg = style.makeLegend(0.50, 0.65, 0.68, 0.8)
        leg.SetTextSize(22)
        cv.SetBottomMargin(0.15)
        cv.SetLeftMargin(0.15)
        cv.SetTopMargin(0.05)
        cv.SetRightMargin(0.05)
        upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
        lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)
        upperPad.SetBottomMargin(0.00001)
        upperPad.SetBorderMode(0)
        upperPad.SetTopMargin(0.05)
        upperPad.SetRightMargin(0.05)
        upperPad.SetLeftMargin(0.13)
        lowerPad.SetTopMargin(0.00001)
        lowerPad.SetBottomMargin(0.4)
        lowerPad.SetRightMargin(0.05)
        lowerPad.SetLeftMargin(0.13)
        lowerPad.SetBorderMode(0)
        #canvas.SetBottomMargin(0.2)
        #canvas.SetTopMargin(0.05)
        upperPad.Draw()
        lowerPad.Draw()
        upperPad.cd()
        model=(f"nominal_{key}", f"nominal_{key}", 20, 0., 1.)
        model_up=(f"up_{key}", f"up_{key}", 20, 0., 1.)
        model_down=(f"down_{key}", f"down_{key}", 20, 0., 1.)

        if uncertainty == "shape":
            rdf = rdf.Define(f"discriminant_{key}_up", f"max(hnlJet_{key}Up_llpdnnx_ratio_LLP_QMU, hnlJet_{key}Up_llpdnnx_ratio_LLP_QE)")
            rdf = rdf.Define(f"discriminant_{key}_down", f"max(hnlJet_{key}Down_llpdnnx_ratio_LLP_QMU, hnlJet_{key}Down_llpdnnx_ratio_LLP_QE)")

            hist_nominal = rdf.Histo1D(model, "discriminant")
            hist_up = rdf.Histo1D(model_up, f"discriminant_{key}_up")
            hist_down = rdf.Histo1D(model_down, f"discriminant_{key}_down")
        else:
            rdf = rdf.Define(f"nominal_{key}", uncertainty)
            rdf = rdf.Define(f"up_{key}",  uncertainty.replace("nominal", "up"))
            rdf = rdf.Define(f"down_{key}", uncertainty.replace("nominal", "down"))


            hist_nominal = rdf.Histo1D(model, "discriminant", f"nominal_{key}")
            hist_up = rdf.Histo1D(model_up, "discriminant", f"up_{key}")
            hist_down = rdf.Histo1D(model_down, "discriminant", f"down_{key}")
        print(uncertainty, hist_up.Integral(), hist_down.Integral(), hist_nominal.Integral())
        difference_percent_up = (hist_up.Integral()-hist_nominal.Integral())/hist_nominal.Integral()*100
        difference_percent_down = (hist_down.Integral()-hist_nominal.Integral())/hist_nominal.Integral()*100


        hist_up.Scale(1./hist_nominal.Integral())
        hist_down.Scale(1./hist_nominal.Integral())
        hist_nominal.Scale(1./hist_nominal.Integral())
        hist_nominal.GetYaxis().SetTitle("Entries")

        hist_nominal.SetMaximum(1.6*max(hist_nominal.GetMaximum(), hist_up.GetMaximum(), hist_down.GetMaximum()))
        hist_nominal.SetLineWidth(2)
        hist_up.SetLineWidth(3)
        hist_down.SetLineWidth(3)
        hist_up.SetLineStyle(2)
        hist_down.SetLineStyle(3)
        hist_up.SetLineColor(color1)
        hist_down.SetLineColor(color2)

        hist_nominal.Draw("HIST")
        hist_up.Draw("HIST SAME")
        hist_down.Draw("HIST SAME")


        leg.AddEntry(f"nominal_{key}", f"{uncertainty_text[key]} (nominal)", "l")
        leg.AddEntry(f"up_{key}", f"{uncertainty_text[key]} (up, {difference_percent_up:.1f} %)", "l")
        leg.AddEntry(f"down_{key}", f"{uncertainty_text[key]} (down, {difference_percent_down:.1f} %)", "l")

        lowerPad.cd()
        axis = hist_nominal.Clone("axis")
        axis.SetMinimum(0.78)
        axis.SetMaximum(1.24)
        if final_state == "mu":
            axis.GetXaxis().SetTitle("P_{q#mu}")
        elif final_state == "e":
            axis.GetXaxis().SetTitle("P_{qe}")

        axis.GetYaxis().SetTitle("Ratio")
        axis.GetYaxis().SetTitleOffset(1.2)
        axis.GetXaxis().SetTitleOffset(3.5)
        axis.Draw("AXIS")

        ratio_up=hist_up.Clone("ratio_up")
        ratio_up.Divide(hist_nominal.Clone())
        ratio_down=hist_down.Clone("ratio_down")
        ratio_down.Divide(hist_nominal.Clone())
        ratio_up.SetLineStyle(1)
        ratio_down.SetLineStyle(1)
        ratio_up.SetMarkerSize(0)
        ratio_down.SetMarkerSize(0)
        ratio_up.Draw("SAME HIST C")
        ratio_down.Draw("SAME HIST C")
        line = ROOT.TLine(0, 1, 1, 1)
        line.Draw("SAME")
        cv.cd()
        leg.Draw("SAME")

        #style.makeCMSText(0.17, 0.93, additionalText="Work in Progress")
        style.makeText(0.15, 0.88, 0.5, 0.88, "m_{N}=10 GeV, c#tau_{0}=1 mm, V_{#mu}=V_{e}=V_{#tau}", size=25)

        cv.SaveAs(f"plots/{final_state}_{key}.pdf")
        cv.SaveAs(f"plots/{final_state}_{key}.png")
