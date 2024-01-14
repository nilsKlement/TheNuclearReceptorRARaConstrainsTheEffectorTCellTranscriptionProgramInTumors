class RNASequenceGeneral:

    # a class describing a default RNA sequence, the actual values given are not specifiable

    def __init__(self, gene, values):
        self.gene = gene
        self.values = values


class RNASequence:

    # a class describing a default RNA sequence, the actual values are known and thus specified

    def __init__(self, gene, fold_change, p_value):
        self.gene = gene
        self.fold_change = fold_change
        self.p_value = p_value