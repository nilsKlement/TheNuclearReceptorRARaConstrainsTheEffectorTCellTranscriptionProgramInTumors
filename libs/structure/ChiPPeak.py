class ChiPPeak:

    # a class representing the peaks contained in a ChiP sequencing

    def __init__(self, gene, peak_start, peak_end, chrom, fold_enrichment, MAXOVERLAPP=5000):
        self.MAXOVERLAPP = MAXOVERLAPP
        self.gene = gene
        self.chrom, self.peak_start, self.peak_end = chrom, int(peak_start), int(peak_end)
        self.fold_enrichment = float(fold_enrichment)

    # defining equal condition to determine if to peaks are overlapping
    def __eq__(self, other):
        if self.chrom == other.chrom:
            if (abs(self.peak_start - other.peak_start) < self.MAXOVERLAPP
                    or abs(self.peak_end - other.peak_end) < self.MAXOVERLAPP):
                return True
        return False
