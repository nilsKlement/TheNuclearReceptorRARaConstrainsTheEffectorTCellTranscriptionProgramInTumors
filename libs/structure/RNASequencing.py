import os, sys
sys.path.insert(1, r"C:\Users\nilsk\OneDrive\Desktop\workspace\publications\RARa_regulated_antitumor_immunity_github\libs\structure")
from RNASequence import RNASequence, RNASequenceGeneral


class RNASequencingGeneral:
    def __init__(self, filename):
        self.name = filename.split("/")[-1][:-4]
        self.sequences, self.sequences_by_name = self.get_sequences(filename)
        self.number_sequences = len(self.sequences)
        self.header = open(filename, "r").readlines()[0][:-2]

    @staticmethod
    def get_sequences(filename):
        raw_data = open(filename, "r").readlines()[1:]
        sequences = []
        sequences_by_name = {}
        for line in raw_data:
            line = line.split(",")
            line[-1] = line[-1][:-2]
            sequences.append(RNASequenceGeneral(line[0], line[1:]))
            sequences_by_name[line[0]] = RNASequenceGeneral(line[0], line[1:])
        return sequences, sequences_by_name


class RNASequencing:
    def __init__(self, filename):
        self.name = filename.split("/")[-1][:-4]
        self.sequences, self.sequences_by_name = self.get_sequences(filename)
        self.number_sequences = len(self.sequences)

    @staticmethod
    def get_sequences(filename):
        raw_data = open(filename, "r").readlines()[1:]
        sequences = []
        sequences_by_name = {}
        for line in raw_data:
            line = line.split(",")
            if line[-1]!="\n" and line[-2]!='#DIV/0!' and line[-1]!='#DIV/0!\n' and line[-2]!='#NUM!' and line[-1]!='#NUM!':
                sequences.append(RNASequence(line[-3], float(line[-2]), float(line[-1])))
                sequences_by_name[line[-3]] = RNASequence(line[-3], float(line[-2]), float(line[-1]))
        return sequences, sequences_by_name

'''
    def divide_by(self, other):
        result = {}
        for sequence in self.sequences:
            if sequence.gen in list(other.sequences_by_name.keys()):
                if other.sequences_by_name[sequence.gen].fold_change != 0:
                    result[sequence.gen] = sequence.fold_change/other.sequences_by_name[sequence.gen].fold_change
        return result

    def get_folds(self, func=lambda x:x):
        result = {}
        for sequence in self.sequences:
            result[sequence.gen] = func(sequence.fold_change)
        return result
'''