import ROOT
import os
from histo import style

color1 = ROOT.kOrange+7
color2 = ROOT.kSpring+3

benchmark_model = "HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03-2017"
ntuple_path = "/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/26Aug21/2017/" 

file_list = ROOT.std.vector('string')()

for f in os.listdir(os.path.join(ntuple_path, benchmark_model)):
    file_list.push_back(os.path.join(ntuple_path, benchmark_model, f))

for final_state in ["qmu", "qe", "q"]:
    rdf = ROOT.RDataFrame("Friends", file_list)
    rdf = rdf.Filter("nominal_dR_l2j<1.3")
    rdf = rdf.Filter("nominal_met<100 and nominal_m_llj>70 and nominal_m_llj<90 and nselectedJets_nominal<5")
    if final_state == "qmu":
        rdf = rdf.Filter("subleadingLeptons_isMuon[0]").Filter("nominal_dR_l2j<0.4")
        rdf = rdf.Define("discriminant", "hnlJet_nominal_llpdnnx_ratio_LLP_QMU")
    elif final_state == "qe":
        rdf = rdf.Filter("subleadingLeptons_isElectron[0]").Filter("nominal_dR_l2j<0.4")
        rdf = rdf.Define("discriminant", "hnlJet_nominal_llpdnnx_ratio_LLP_QE")
    elif final_state == "q":
        rdf = rdf.Filter("nominal_dR_l2j>0.4")
        rdf = rdf.Define("discriminant", "hnlJet_nominal_llpdnnx_ratio_LLP_Q")

    rdf = rdf.Define("prefired", "nselectedL1PreFiringJets_nominal>0")

    cv = style.makeCanvas()
    leg = style.makeLegend(0.50, 0.65, 0.68, 0.8)
    leg.SetTextSize(22)


    hist_all = rdf.Histo1D(("tagger", "tagger", 10, 0.2, 1), "discriminant")
    hist_prefired = rdf.Histo1D(("tagger", "tagger", 10, 0.2, 1), "discriminant", "prefired")

    hist_all=hist_all.Clone()
    hist_prefired=hist_prefired.Clone()

    cv.Draw("")
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.15)

    efficiency = ROOT.TEfficiency(hist_prefired, hist_all)
    efficiency.Draw("")

    if final_state == "qmu":
        efficiency.SetTitle(";P(q#mu);fraction prefired")
    elif final_state == "qe":
        efficiency.SetTitle(";P(qe);fraction prefired")
    elif final_state == "q":
        efficiency.SetTitle(";P(q);fraction prefired")
    style.makeText(0.15, 0.88, 0.5, 0.88, "m_{N}=10 GeV, c#tau_{0}=1 mm, V_{#mu}=V_{e}=V_{#tau}", size=25)

    cv.SaveAs(f"{final_state}.pdf")
    cv.SaveAs(f"{final_state}.png")
