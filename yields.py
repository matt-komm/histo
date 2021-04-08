import pandas as pd
import numpy as np
import json
import style
import ROOT
import math
from scipy.stats import chi2

# Minimal "ABCD" estimator
#
# m.mieskolainen@imperial.ac.uk, 2021
def ABCD_eq(a,b,c):
    """
    Basic formula
    """
    return a*c/b


def ABCD_err(a,b,c):
    """ 
    Analytical error propagation (1st order Taylor expansion)
    """
    return np.sqrt( (c/b)**2 * a + (a*c/b**2)**2 * b + (a/b)**2 * c)
    #return ABCD_eq(a,b,c) * np.sqrt(1/a + 1/c + 1/b) # Equivalent


def test_abcd():
        
    # ---------------------------
    # INPUT DATA
    A = 100
    B = 50
    C = 1
    # ---------------------------

    # Simple estimate
    D       = ABCD_eq(a=A,b=B,c=C)
    sigma_D = ABCD_err(a=A,b=B,c=C)

    # Random sample ~ Bootstrap via implicit/explicit Poisson assumption
    N     = int(1e6)
    A_new = np.random.poisson(lam=A, size=N)
    B_new = np.random.poisson(lam=B, size=N)
    C_new = np.random.poisson(lam=C, size=N)
    D_new = ABCD_eq(a=A_new, b=B_new, c=C_new)

    # Confidence interval via percentile bootstrap
    alpha = 1 - 0.68 # confidence level

    sigma_D_BS_lo = np.percentile(D_new, 100*(alpha/2))
    sigma_D_BS_hi = np.percentile(D_new, 100*(1-alpha/2))

    print(f'D = {D}, errorprop_CL = [{D-sigma_D:0.2f}, {D+sigma_D:0.2f}], bootstrap_CL = [{sigma_D_BS_lo:0.2f}, {sigma_D_BS_hi:0.2f}]')


def make_pretty(s: str) -> str:
    """ Make pretty the category text"""
    s = s.replace("mumu", "$\mu\mu$")
    s = s.replace("mue", "$\mu e$")
    s = s.replace("emu", "$e\mu$")
    s = s.replace("ee", "$ee$")
    s = s.replace("_", ",")
    s = s.replace("displaced", "displaced")
    s = s.replace("prompt", "prompt")
    return s


def get_hist(file_name: str, hist_name: str) -> ROOT.TH1D:
    rootFile = ROOT.TFile(file_name)
    hist = rootFile.Get(hist_name)
    hist = hist.Clone()
    hist.SetDirectory(0)
    rootFile.Close()
    return hist


def poisson_interval(k:float, alpha:float=0.3173) -> float: 
    """
    uses chisquared info to get the poisson interval. Uses scipy.stats 
    (imports in function). 
    """
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    if k == 0: 
        low = 0.0
    return low, high


def prediction(A: float, B: float, C:float) -> str:

    # Mikael Mieskolainen
    # Random sample ~ Bootstrap via implicit/explicit Poisson assumption
    N = int(1e6)
    A_new = np.random.poisson(lam=A, size=N)
    B_new = np.random.poisson(lam=B, size=N)
    C_new = np.random.poisson(lam=C, size=N)
    D_new = ABCD_eq(a=A_new, b=B_new, c=C_new)

    # Confidence interval via percentile bootstrap
    alpha = 1 - 0.68 # confidence level

    sigma_D_BS_lo = np.percentile(D_new, 100*(alpha/2))
    sigma_D_BS_hi = np.percentile(D_new, 100*(1-alpha/2))
    D_new = np.mean(D_new)

    err_D_lower = D_new - sigma_D_BS_lo
    err_D_upper = sigma_D_BS_hi - D_new

    return D_new, err_D_lower, err_D_upper


def plot_yields(df: pd.DataFrame, topology: str="merged", years: list=["2016", "2017", "2018"]) -> None:

    """
    Make a plot of a pandas dataframe
    """

    for year in years:
        bkg_predictions = df[f"pred{year}"].tolist()
        sig_predictions = df[f"signal{year}"].tolist()
        bin_names = list(df)
        n_bins_total = len(sig_predictions)
        print(sig_predictions, bkg_predictions, bin_names, n_bins_total)

        bkg_hist = ROOT.TH1D(year, year, n_bins_total, -0.5, n_bins_total-0.5)
        sig_hist = bkg_hist.Clone(f"bkg_{year}")

        for j, (bin_name, sig_pred, bkg_pred) in enumerate(zip(bin_names, sig_predictions, bkg_predictions)):
            # if j == 0:
            #     bin_name = replace_name(category_name)
            # if j == 1:
            #     bin_name = ""
            print(j, bin_name, sig_pred, bkg_pred)
            bkg_hist.GetXaxis().SetBinLabel(j, bin_name)
            bkg_hist.GetXaxis().LabelsOption("v")
            bkg_hist.SetBinContent(j, bkg_pred)
            sig_hist.SetBinContent(j, sig_pred)
        cv = style.makeCanvas()
        cv.Draw("")
        bkg_hist.Draw("")
        sig_hist.Draw("SAME")
        cv.SaveAs(f"yields/{topology}_{year}.pdf")


category_file = 'categories_2l.json'

with open(category_file, 'r') as fp:
    category_dict = json.load(fp)


for topology in ["merged", "resolved"]:
    if topology == "merged":
        idx = 1
    else:
        idx = 2
    category_names = []
    yields = {"2016": [], "2017": [], "2018": []}
    signal_yields = {"2016": [], "2017": [], "2018": []}

    for category in category_dict.keys():
        category_names.append(make_pretty(category))
        for year in yields.keys():
            hist_A = get_hist(f"hists_merged/{year}.root", f"{category}_A/data")
            hist_B = get_hist(f"hists_merged/{year}.root", f"{category}_B/data")
            hist_C = get_hist(f"hists_merged/{year}.root", f"{category}_C/data")
            hist_signal = get_hist(f"hists_merged/HNL_majorana_all_ctau1p0e00_massHNL8p0_Vall2p119e-03_{year}.root", f"{category}_D/HNL_coupling_7")
            a = hist_A.GetBinContent(idx)
            b = hist_B.GetBinContent(idx)
            c = hist_C.GetBinContent(idx)
            d, d_up, d_low = prediction(a, b, c)
            d_string =  f'${d:.1f}^{{+{d_up:.1f}}}_{{-{d_low:.1f}}}$'

            yields[year].append(d_string)
            hist_signal.Scale(0.1932398841208851)

            signal_yield = hist_signal.GetBinContent(idx)
            signal_uncertainty = hist_signal.GetBinError(idx)
            signal_string = f"${signal_yield:0.1f} \\pm {signal_uncertainty:0.1f}$"

            signal_yields[year].append(signal_string)

    df = pd.DataFrame(dict(name=category_names, pred2016=yields["2016"], pred2017=yields["2017"], pred2018=yields["2018"], signal2016=signal_yields["2016"], signal2017=signal_yields["2017"], signal2018=signal_yields["2018"]))
    plot_yields(df, topology=topology)
    print(df.to_latex(index=False, escape=False, caption=f"Expected background and signal yields for {topology} categories. The signal model shown is Majorana HNL of $8\\GeV, c\\tau_0=1\\unit{{mm}}, V_\\mu=V_e$.", label=f"tab:yields_{topology}", column_format='c|c|c|c|c|c|c'))
