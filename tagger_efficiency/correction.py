import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import ROOT
from scipy.optimize import leastsq


def read_hist(category, year):
    name = "hnlJets_svMatchedTracks_adapted_nominal_dxy"
    f = ROOT.TFile.Open(f"plots_tagger_eff/{category}_high_met_hnlJets_svMatchedTracks_adapted_nominal_dxy_{year}.root")
    canvas = f.Get(name)    
    upperPad = canvas.GetPrimitive("upperPad")
    stack = upperPad.GetPrimitive(name)
    h_data = upperPad.GetPrimitive("data")
    h_mc = stack.GetStack().Last()
    h_mc.Scale(h_data.Integral()/h_mc.Integral())
    h_ratio = h_data.Clone("ratio")
    h_ratio.Divide(h_mc)

    centers = []
    ratios = []
    errors = []
    for i in range(h_ratio.GetNbinsX()):
        center = h_mc.GetBinCenter(i+1)
        ratio = h_ratio.GetBinContent(i+1)
        error = h_ratio.GetBinError(i+1)
        centers.append(center)
        ratios.append(ratio)
        errors.append(error)
    f.Close()
    return centers, ratios, errors

for year in ["2016", "2017", "2018"]:

    centers_mu, ratios_mu, errors_mu = read_hist("mu1", year)
    centers_e, ratios_e, errors_e = read_hist("el1",year)

    centers = np.asarray(centers_mu + centers_e)
    ratios = np.asarray(ratios_mu + ratios_e)
    errors = np.asarray(errors_mu + errors_e)

    x = np.log10(np.geomspace(centers_mu[0], centers_mu[-1]))
    log_centers = np.log10(centers)
    def fitfunc(p, x):
        y = np.piecewise(x, [x < p[0], x >= p[0]],
                        [lambda x:p[2]*x + p[1]-p[2]*p[0], lambda x:p[3]*x + p[1]-p[3]*p[0]])
        return y
        
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    pinit = [-1.5,  1,  0.1, -0.1]

    N = int(1e4)
    n = len(centers)
    bootstrap_indices = np.random.randint(0, n, size=(N, n))
    bootstrap_parameters = np.zeros((N, len(pinit)))

    for i, indices in enumerate(bootstrap_indices):
        bootstrap_centers = log_centers[indices]
        bootstrap_ratios = ratios[indices]
        bootstrap_errors = errors[indices]
        pvals = leastsq(errfunc, pinit, args=(bootstrap_centers, bootstrap_ratios, bootstrap_errors))
        bootstrap_parameters[i] = pvals[0]

    pfinal = np.mean(bootstrap_parameters, axis=0)
    pstd = np.std(bootstrap_parameters, axis=0)
    cov = np.cov(bootstrap_parameters.T)

    #print(pfinal, pstd)
    #print(cov)

    preds = fitfunc(pfinal, log_centers)
    residuals = np.divide((preds - ratios) ** 2, errors)
    chi_squared = np.sum(residuals)

    plt.style.use(hep.style.ROOT)
    fig, ax = plt.subplots()
    plt.errorbar(centers_mu, ratios_mu, yerr=errors_mu, label='points (muon)')
    plt.errorbar(centers_e, ratios_e, yerr=errors_e, label='points (electron)')
    plt.plot(np.power(10, x), fitfunc(pfinal, x), label=r"2nd-order fit, $\chi^2=$"+f"{chi_squared:.2f}")
    ax.legend()
    ax.set_xscale('log')
    ax.set_xlabel(r'track $d_{xy}$ [cm]')
    ax.set_ylabel('Correction (Data/MC)')
    plt.savefig(f'tagger_efficiency/fits/combined_{year}.pdf')
