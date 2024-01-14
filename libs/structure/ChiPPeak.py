import numpy as np

MAXOVERLAPP = 5000

class ChiPPeak:
    def __init__(self, gene, peak_start, peak_end, chrom, fold_enrichment):
        self.gene = gene
        self.chrom, self.peak_start, self.peak_end = chrom, int(peak_start), int(peak_end)
        self.fold_enrichment = float(fold_enrichment)

    def __str__(self):
        return " ".join([str(element) for element in [self.gene, self.chrom, self.peak_start, self.peak_end, self.fold_enrichment]])

    def __eq__(self, other):
        if self.chrom == other.chrom:
            if (abs(self.peak_start - other.peak_start) < MAXOVERLAPP \
                    or abs(self.peak_end - other.peak_end) < MAXOVERLAPP):
                return True
        return False