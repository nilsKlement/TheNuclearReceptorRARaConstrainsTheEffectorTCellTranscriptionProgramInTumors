class RNASequenceGeneral:
    def __init__(self, gene, values):
        self.gene = gene
        self.values = values

class RNASequence:
    def __init__(self, gene, fold_change, p_value):
        self.gene = gene
        self.fold_change = fold_change
        self.p_value = p_value