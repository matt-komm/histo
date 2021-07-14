import uproot
import pandas as pd
from pathlib import Path
import ROOT
from histo import style
style.makeColorTable()

frames = []
ntuple_path = Path("/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/29Jun21/2016/HNL_dirac_pt20_ctau1p0e00_massHNL10p0_Vall1p664e-03-2016")

for f in ntuple_path.iterdir():
    events = uproot.open(f)["Friends"]
    frames.append(events.pandas.df(["hnlJet_nominal_isTrue*", "nominal_dR_l2j"]))
df = pd.concat(frames)
df = df.query('nominal_dR_l2j<1.3 ')
df = df.query('hnlJet_nominal_isTrueQ>0 or hnlJet_nominal_isTrueQE>0 or hnlJet_nominal_isTrueQTAU>0 or hnlJet_nominal_isTrueQMU>0')
df['true_merged'] = df['hnlJet_nominal_isTrueQE']+df['hnlJet_nominal_isTrueQMU']
df['true_resolved'] = df['hnlJet_nominal_isTrueQ']+df['hnlJet_nominal_isTrueQTAU']
df['reco_merged'] = df['nominal_dR_l2j']<0.4
df['reco_resolved'] = df['nominal_dR_l2j']>0.4

merged_merged = len(df.query('reco_merged>0 and true_merged>0'))/len(df.query('true_merged>0'))*100
merged_resolved = len(df.query('reco_resolved>0 and true_merged>0'))/len(df.query('true_merged>0'))*100
resolved_merged = len(df.query('reco_merged>0 and true_resolved>0'))/len(df.query('true_resolved>0'))*100
resolved_resolved = len(df.query('reco_resolved>0 and true_resolved>0'))/len(df.query('true_resolved>0'))*100

print(merged_merged, merged_resolved, resolved_merged, resolved_resolved)

hist = ROOT.TH2D("reco_gen", "reco_gen", 2, -0.5, 1.5, 2, -0.5, 1.5)
hist.Fill(0, 0, merged_merged)
hist.Fill(1, 0, merged_resolved)
hist.Fill(0, 1, resolved_merged)
hist.Fill(1, 1, resolved_resolved)
hist.GetXaxis().SetBinLabel(1, "merged (reco)")
hist.GetXaxis().SetBinLabel(2, "resolved (reco)")
hist.GetYaxis().SetBinLabel(1, "LLP_QL")
hist.GetYaxis().SetBinLabel(2, "LLP_Q")
hist.GetZaxis().SetTitle("Fraction (%)")
hist.SetMarkerSize(3)
cv = style.makeCanvas()
cv.SetBottomMargin(0.12)
cv.SetLeftMargin(0.15)
cv.SetRightMargin(0.18)
cv.Draw("")
hist.Draw("COLZ text")
cv.SaveAs("reco_gen.pdf")
cv.SaveAs("reco_gen.png")