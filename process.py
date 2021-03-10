import ROOT
from . import style
import json

with open("/vols/cms/LLP/color_dict_mk.json") as json_file:
    color_dict = json.load(json_file)

# A process is a combination of several "Samples" which are all added up internally
class Process:
    def __init__(self, name, title, linecolor=None, fillcolor=None):
        self.name = name
        self.title = title
        self.linecolor = linecolor
        self.fillcolor = fillcolor
        if not linecolor:
            self.linecolor = color_dict[name][1]
        if not fillcolor:
            self.fillcolor = color_dict[name][0]
        self.hists = []
        self.rdfs = []

    def add(self, *args):
        for arg in args:
            self.rdfs.append(arg.rdf)

    def Define(self, var, varexp):
        for i, rdf in enumerate(self.rdfs):
            if var not in rdf.GetColumnNames():
                self.rdfs[i] = rdf.Define(var, varexp)

    def Histo1D(self, args, varexp, weight="weightNominal", cut=None):
        for i, rdf in enumerate(self.rdfs):
            if cut is not None:
                _rdf = rdf.Filter(cut)
            else:
                _rdf = rdf
            if i == 0:
                hist = _rdf.Histo1D(args, varexp, weight)
                hist = hist.Clone()
                hist.Sumw2()
            else:
                tmp_hist = _rdf.Histo1D(args, varexp, weight)
                tmp_hist = tmp_hist.Clone()
                tmp_hist.Sumw2()
                hist.Add(tmp_hist)

        hist.GetXaxis().SetTitle(args[1])
        hist.SetLineColor(ROOT.TColor.GetColor(self.linecolor))
        hist.SetFillColor(ROOT.TColor.GetColor(self.fillcolor))
        hist.SetLineWidth(2)
        hist.SetName(self.name)
        hist.SetTitle(self.title)
        self.hists.append(hist)
        return self.hists[-1]
    
    def AsNumpy(self, *args, cut=None):
        output_dict = {}
        for _, rdf in enumerate(self.rdfs):
            _rdf = rdf.Filter(cut)
            _dict = _rdf.AsNumpy(columns=args)
            output_dict = {**output_dict, **_dict}
        return output_dict

    def Histo2D(self, args, varexp1, varexp2, weight="weight", cut=None):
        for i, rdf in enumerate(self.rdfs):
            if cut is not None:
                _rdf = rdf.Filter(cut)
            else:
                _rdf = rdf
            if i == 0:
                hist = _rdf.Histo2D(args, varexp1, varexp2, weight)
                hist = hist.Clone()
                hist.Sumw2()
            else:
                tmp_hist = _rdf.Histo2D(args, varexp1, varexp2, weight)
                tmp_hist = tmp_hist.Clone()
                tmp_hist.Sumw2()
                hist.Add(tmp_hist)
        self.hists.append(hist)
        return self.hists[-1]
