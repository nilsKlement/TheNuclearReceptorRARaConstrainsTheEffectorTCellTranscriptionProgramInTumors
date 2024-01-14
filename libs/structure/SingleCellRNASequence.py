class SingleCellRNASequenceGeneral:
    def __init__(self, gene, values):
        self.gene = gene
        self.values = values


class SingleCellRNASequence:
    def __init__(self, gene, WT_Log2FoldChange, KO_Log2FoldChange, TG_Log2FoldChange):
        self.gene = gene
        self.WT_Log2FoldChange = WT_Log2FoldChange
        self.KO_Log2FoldChange = KO_Log2FoldChange
        self.TG_Log2FoldChange = TG_Log2FoldChange