import pandas as pd
import os, sys
sys.path.insert(1, r"C:\Users\nilsk\OneDrive\Desktop\workspace\publications\RARa_regulated_antitumor_immunity\libs\structure")
from ChipGen import ChiPGen


class ChiPSeq:
    GEN = "symbol"
    CHROMBROWSER = "chromBrowser"
    FOLDENRICHMENT = "foldEnrichment"
    MAXOVERLAPP = 5000

    def __init__(self, files: list, name=None, no_latex_name=None, color=None):
        self.name, self.color = name, color
        if no_latex_name is None:
            self.no_latex_name = name
        else:
            self.no_latex_name = no_latex_name
        self.gens = []
        for file in files:
            self.add_measurement(file)
        for gen in self.gens:
            gen.average_fold_enrichments()
        self.number_gens = len(self.gens)

        self.save_to_csv()

    def add_measurement(self, file):
        df = pd.read_csv(file)
        gens = df.loc[:, self.GEN].to_list()
        bindings = df.loc[:, self.CHROMBROWSER].to_list()
        fold_enrichments = df.loc[:, self.FOLDENRICHMENT].to_list()
        for gen, binding, fold_enrichment in zip(gens, bindings, fold_enrichments):
            temp_gen = ChiPGen(gen, binding, fold_enrichment)
            if not isinstance(gen, float):
                if self.gen_unique(temp_gen):
                    self.gens.append(temp_gen)

    def gen_unique(self, gen):
        for _gen in self.gens:
            if _gen == gen:
                _gen.add_fold_enrichment(gen.fold_enrichments[0])
                return False
        return True

    def __str__(self):
        return "Chip seq with "+str(self.number_gens)+" sequenced gens"

    def save_to_csv(self):
        table = []
        for gen in self.gens:
            table.append([gen.gene, gen.fold_enrichment, gen.chrom, gen.peak_start, gen.peak_end])
        pd.DataFrame(table).to_csv("results/data/combainedChiP___"+self.no_latex_name+ ".csv", index=False, header=False)