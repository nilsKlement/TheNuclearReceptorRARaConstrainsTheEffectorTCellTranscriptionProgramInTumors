import pandas as pd
import os, sys
sys.path.insert(1, r"C:\Users\nilsk\OneDrive\Desktop\workspace\publications\RARa_regulated_antitumor_immunity_github\libs\structure")
from ChiPPeak import ChiPPeak


class ChiPSeq:

    # a class representing a ChiP sequencing

    # defining the headers of the measurement files
    GEN = "Symbol"
    CHROMBROWSER1 = "start"
    CHROMBROWSER2 = "end"
    CHROMOSOM = "chr"
    FOLDENRICHMENT = "fold_enrichment"
    MAXOVERLAPP = 5000

    def __init__(self, files: list, name=None, no_latex_name=None, color=None):
        # variables for a proper visualization
        self.name, self.color = name, color
        if no_latex_name is None:
            self.no_latex_name = name
        else:
            self.no_latex_name = no_latex_name

        self.peaks = []
        # iteration over several files, it measurement is divided in multiple files
        for file in files:
            self.add_measurement(file)

        self.number_gens = len(self.peaks)

    def add_measurement(self, file):
        # loading the measurement file and dividing into lists for every characteristic property
        df = pd.read_csv(file)
        genes = df.loc[:, self.GEN].to_list()
        starts = df.loc[:, self.CHROMBROWSER1].to_list()
        ends = df.loc[:, self.CHROMBROWSER2].to_list()
        chromosoms = df.loc[:, self.CHROMOSOM].to_list()
        fold_enrichments = df.loc[:, self.FOLDENRICHMENT].to_list()

        # create an instance of ChipPeak for each peak in the file
        # by iterating through them based on their defined characteristic properties
        for gen, start, end, chromosom, fold_enrichment in zip(genes, starts, ends, chromosoms, fold_enrichments):

            # empty cells are treated as nan's initially, replacing these by an empty string
            if isinstance(gen, float):
                gen = ""

            self.peaks.append(ChiPPeak(gen, start, end, chromosom, fold_enrichment, MAXOVERLAPP=self.MAXOVERLAPP))
