import pandas as pd
import os, sys
sys.path.insert(1, r"C:\Users\nilsk\OneDrive\Desktop\workspace\publications\RARa_regulated_antitumor_immunity_github\libs\structure")
from ChiPPeak import ChiPPeak
import numpy as np


class ChiPSeq:
    GEN = "Symbol"
    CHROMBROWSER1 = "start"
    CHROMBROWSER2 = "end"
    CHROMOSOM = "chr"
    FOLDENRICHMENT = "fold_enrichment"
    MAXOVERLAPP = 5000

    def __init__(self, files: list, name=None, no_latex_name=None, color=None):
        self.name, self.color = name, color
        if no_latex_name is None:
            self.no_latex_name = name
        else:
            self.no_latex_name = no_latex_name
        self.peaks = []
        for file in files:
            self.add_measurement(file)
        self.number_gens = len(self.peaks)

    def add_measurement(self, file):
        df = pd.read_csv(file)
        genes = df.loc[:, self.GEN].to_list()
        starts = df.loc[:, self.CHROMBROWSER1].to_list()
        ends = df.loc[:, self.CHROMBROWSER2].to_list()
        chromosoms = df.loc[:, self.CHROMOSOM].to_list()
        fold_enrichments = df.loc[:, self.FOLDENRICHMENT].to_list()
        for gen, start, end, chromosom, fold_enrichment in zip(genes, starts, ends, chromosoms, fold_enrichments):
            if isinstance(gen, float):
                gen = ""
            self.peaks.append(ChiPPeak(gen, start, end, chromosom, fold_enrichment))