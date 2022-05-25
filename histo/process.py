import ROOT
import numpy as np
import pandas as pd

color_dict = {
    "wjets": ["#54bf59", "#327235"],
    "dyjets": ["#4cbff2", "#2d7291"],
    "vgamma": ["#ba3ff4", "#6f2692"],
    "topbkg": ["#ffcc00", "#997a00"],
    "qcd": ["#d8d8d8", "#828282"],
    "data": ["#000000","#000000"]
}

class Process:
    """ 
    The process class is a combination of several "Samples" (with different cross-section etc) 
    which are all added up internally for plotting/histogramming
    Each sample is represented as a ROOT RDataFrame
    """
    def __init__(self, name:str, title:str, linecolor:str=None, fillcolor:str=None):
        self.name = name
        self.title = title
        self.linecolor = linecolor
        self.fillcolor = fillcolor
        if not linecolor:
            if name in color_dict:
                self.linecolor = color_dict[name][1]
            else:
                self.linecolor = '#ffffff'
        if not fillcolor:
            if name in color_dict:
                self.fillcolor = color_dict[name][0]
            else:
                self.fillcolor = '#ffffff'
        
        self.hists = []
        self.rdfs = []

    def Add(self, *args):
        """ Add sample(s)"""
        for arg in args:
            self.rdfs.append(arg.rdf)

    def Define(self, var:str, varexp:str):
        """ Define a variable with name var and expression varexp"""
        if len(self.rdfs) == 0:
            pass
        if var not in self.rdfs[0].GetColumnNames():
            self.rdfs = [rdf.Define(var, varexp) for rdf in self.rdfs]

    def Histo1D(self, args:list, varexp:str, weight:str="weightNominal", cut:str=None):
        """ 
        Produce a 1D histogram
        Takes: 
        args: pass to RDataFrame (bins etc..)
        varexp: variable to plot
        weight: weight to apply
        cut: cut to apply
        Returns: ROOT TH1D Histogram

        """
        for i, rdf in enumerate(self.rdfs):
            if cut is not None:
                _rdf = rdf.Filter(cut)
            else:
                _rdf = rdf
            if i == 0:
                hist = _rdf.Histo1D(args, varexp, weight)
                hist = hist.Clone()
                hist.SetDirectory(0)
                hist.Sumw2()
            else:
                tmp_hist = _rdf.Histo1D(args, varexp, weight)
                tmp_hist = tmp_hist.Clone()
                tmp_hist.SetDirectory(0)
                tmp_hist.Sumw2()
                hist.Add(tmp_hist)

        hist.GetXaxis().SetTitle(args[1])
        hist.SetLineColor(ROOT.TColor.GetColor(self.linecolor))
        hist.SetFillColor(ROOT.TColor.GetColor(self.fillcolor))
        hist.SetLineWidth(2)
        hist.SetName(self.name)
        hist.SetTitle(self.title)
        self.hists.append(hist) # To ensure no memory leaks
        return self.hists[-1]
 
    def Histo2D(self, args: list, varexp1: str, varexp2: str, weight: str="weight", cut: str=None):
        """ 
        Produce a 2D histogram
        Takes: 
        args: pass to RDataFrame (bins etc..)
        varexp1: variable x to plot
        varexp2: variable x to plot
        weight: weight to apply
        cut: cut to apply
        Returns: ROOT TH1D Histogram

        """
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
    
    def AsNumpy(self, *args):
        """ Return the stored dataframe as numpy array """
        output = []
        for rdf in self.rdfs:
            out = rdf.AsNumpy(columns=args)
            out = pd.DataFrame(out, columns=args)
            output.append(out)

        return pd.concat(output)
