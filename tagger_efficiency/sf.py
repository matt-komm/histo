import pickle as pkl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import sys
import os
sys.path.append('tagger_efficiency/covidgen/covidgen')
from estimators import clopper_pearson_err

class Cut():
    def __init__(self, cutstring, label, title):
        self.cutstring = cutstring
        self.label = label
        self.title = title

def weight_array(r):
    return  np.ones(r.length) * r.weight

def histogram2d(x, y, bins_x, bins_y, **kwargs):
    hist, bins_pt, bins_eta = np.histogram2d(x, y, bins=(bins_x, bins_y), **kwargs)
    return hist

def plot1d(arrays, labels, bins, title, x_label, log_x):
    plt.style.use(hep.style.ROOT)
    fig, ax = plt.subplots()
    if log_x:
        ax.set_xscale('log')
    #hep.cms.label(loc=0)
    for array, label in zip(arrays, labels):
        h, bins = np.histogram(array, bins, density=True)
        print(h, bins)
        hep.histplot(h, bins, label=label)
    ax.set_xlabel(x_label)
    ax.legend()
    plt.savefig(f"tagger_efficiency/{title}.pdf")

df_data = pd.concat([pd.read_pickle(r'tagger_efficiency/output_20Apr21/2016/'+x+".pkl") for x in ["muon", "electron"]])
df_bkg = pd.concat([pd.read_pickle(r'tagger_efficiency/output_20Apr21/2016/'+x+".pkl") for x in ["dyjets", "topbkg", "wjets", "vgamma"]])
df_bkg['length'] = df_bkg['track_pt'].map(lambda x: len(x))
df_bkg['weight'] = df_bkg.apply(weight_array, axis=1)
print(df_bkg)


merged_mu = Cut("merged and muon", "merged_mu", r"merged, $\mu$")
merged_e = Cut("merged and electron", "merged_e", "merged, e")
resolved_mu = Cut("resolved and muon", "resolved_mu", r"resolved, $\mu$")
resolved_e = Cut("resolved and electron", "resolved_e", "resolved, e")

#labels = [r"DY+jets, $m_{\ell\ell} > 80 GeV$", r"data, $m_{\ell\ell} > 80 GeV$", r"HNL, $m_{N} = 10 GeV, c\tau_0 = 1 mm$" ]
#plot1d([df_bkg['track_pt'], df_data['track_pt'], df_sig['track_pt']], labels, np.geomspace(1., 20., num=25), "pt", x_label=r"leading track $p_\mathrm{T}$ [GeV]", log_x=True)
#plot1d([df_bkg['track_dxy'], df_data['track_dxy'], df_sig['track_dxy']], labels, np.geomspace(0.01, 0.3, num=25), "dxy", x_label=r"leading track $d_{xy}$ [cm]", log_x=True)
for cut in [merged_mu, merged_e, resolved_mu, resolved_e]:
    plt.style.use(hep.style.ROOT)

    bins_pt = [0.5, 1, 3, 8, 20]
    bins_dxy = np.geomspace(0.01, 20., num=5)

    track_pt = np.concatenate(df_bkg['track_pt'].to_numpy().flatten())
    track_dxy = np.concatenate(df_bkg['track_dxy'].to_numpy().flatten())
    weights = np.concatenate(df_bkg['weight'].to_numpy().flatten())

    track_pt_data = np.concatenate(df_data['track_pt'].to_numpy().flatten())
    track_dxy_data = np.concatenate(df_data['track_dxy'].to_numpy().flatten())

    hist = histogram2d(track_pt, track_dxy, bins_pt, bins_dxy, weights=weights)
    hist_data = histogram2d(track_pt_data, track_dxy_data, bins_pt, bins_dxy)

    fig, ax = plt.subplots()
    fig.suptitle(cut.title, fontsize=30)
    hep.cms.label(loc=0)
    hep.hist2dplot(hist_data, bins_pt, bins_dxy, labels=True)

    ax.set_xticklabels(labels=bins_pt, rotation=45)
    ax.set_xlabel(r'leading track $p_\mathrm{T}$ [GeV]')
    ax.set_ylabel(r'leading track $|d_{xy}| [cm]$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.tight_layout()
    plt.show()
    plt.savefig(f'tagger_efficiency/counts_numerator_data_{cut.label}.pdf')

    hist_ratio = []
    labels = []
    for nums, denoms in zip(hist_data, hist):
        row = []
        label_row = []
        for num, denom in zip(nums, denoms):
            eff = np.divide(num, denom)
            sf = num/denom
            print(sf)
            #print(sf, sf_CI)
            #sf_err_down = sf-sf_CI[0]
            #sf_err_up = sf_CI[1]-sf
            row.append(sf)
        hist_ratio.append(row)
        #label_row.append(f"${sf:0.2f}^{{{+sf_err_up:0.2f}}}_{{{-sf_err_down:0.2f}}}$")
        #labels.append(label_row)

    fig, ax = plt.subplots()
    fig.suptitle(cut.title, fontsize=30)
    hep.cms.label(loc=0)
    labels = np.asarray(labels)
    hep.hist2dplot(hist_ratio, bins_pt, bins_dxy)

    _xbin_centers = bins_pt[1:] - np.diff(bins_pt) / float(2)
    _ybin_centers = bins_dxy[1:] - np.diff(bins_dxy) / float(2)
    for ix, xc in enumerate(_xbin_centers):
        for iy, yc in enumerate(_ybin_centers):
            color = "black"
            ax.text(xc, yc, labels[ix, iy], ha="center", va="center", color=color, size=10)

    ax.set_xlabel(r'leading track $p_\mathrm{T}$ [GeV]')
    ax.set_xticklabels(labels=bins_pt, rotation=45)
    ax.set_ylabel(r'leading track $|d_{xy}| [cm]$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()
    plt.savefig(f'tagger_efficiency/{cut.label}.pdf')
