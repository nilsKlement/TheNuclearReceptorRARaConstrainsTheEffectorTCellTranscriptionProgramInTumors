import numpy as np


MAXOVERLAPP = 5000


class ChiPGen:
    def __init__(self, gen, chrom_browser, fold_enrichment):
        self.gen = gen
        self.fold_enrichments = [fold_enrichment]
        self.chrom, self.binding_start, self.binding_end = self.encode_chrombrowser(chrom_browser)
        self.fold_enrichment, self.fold_enrichment_std = None, None

    @staticmethod
    def encode_chrombrowser(chrom_browser):
        chromosom = chrom_browser.split(":")[0]
        int1, int2 = chrom_browser.split(":")[-1].split("-")
        return chromosom, int(int1), int(int2)

    def add_fold_enrichment(self, fold_enrichment):
        self.fold_enrichments.append(fold_enrichment)

    def average_fold_enrichments(self):
        self.fold_enrichment = np.mean(self.fold_enrichments)
        self.fold_enrichment_std = np.std(self.fold_enrichments)

    def __str__(self):
        return " ".join([str(element) for element in [self.gen, self.chrom, self.binding_start, self.binding_end, self.fold_enrichment]])

    def __eq__(self, other):
        if self.chrom == other.chrom:
            if (abs(self.binding_start - other.peak_start) < MAXOVERLAPP \
                    or abs(self.binding_end - other.peak_end) < MAXOVERLAPP):
                return True
        return False