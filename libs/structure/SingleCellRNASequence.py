class SingleCellRNASequenceGeneral:

    # a class describing a default scRNA sequence, the actual values given are not specifiable

    def __init__(self, gene, values):
        self.gene = gene
        self.values = values


class SingleCellRNASequence:

    # a class describing a default scRNA sequence, the actual values are known and thus specified

    def __init__(self, gene, WT_Log2FoldChange, KO_Log2FoldChange, TG_Log2FoldChange):
        self.gene = gene
        self.WT_Log2FoldChange = WT_Log2FoldChange
        self.KO_Log2FoldChange = KO_Log2FoldChange
        self.TG_Log2FoldChange = TG_Log2FoldChange