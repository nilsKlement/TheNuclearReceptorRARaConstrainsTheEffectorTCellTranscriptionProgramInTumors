import os, sys
sys.path.insert(1, r"C:\Users\nilsk\OneDrive\Desktop\workspace\publications\RARa_regulated_antitumor_immunity_github\libs\structure")
from SingleCellRNASequence import SingleCellRNASequence, SingleCellRNASequenceGeneral


class SingleCellRNASequencingGeneral:

    # a class based on the general scRNA sequence for a scRNA sequencing

    def __init__(self, filename):
        self.name = filename.split("/")[-1][:-4]
        self.sequences, self.sequences_by_name = self.get_sequences(filename)
        self.number_sequences = len(self.sequences)
        self.header = open(filename, "r").readlines()[0][:-2]

    # returns sequences as list - for iterations - and as dict - for access by name
    @staticmethod
    def get_sequences(filename):
        # loading the file
        raw_data = open(filename, "r").readlines()[1:]
        sequences = []
        sequences_by_name = {}
        for line in raw_data:
            # delimiter in the file is ",", returns list of properties as strings
            line = line.split(",")
            line[-1] = line[-1][:-2]
            # appending general RNA sequence
            sequences.append(SingleCellRNASequenceGeneral(line[0], line[1:]))
            sequences_by_name[line[0]] = SingleCellRNASequenceGeneral(line[0], line[1:])
        return sequences, sequences_by_name


class SingleCellRNASequencing:

    # a class based on the specific scRNA sequence for a scRNA sequencing

    def __init__(self, filename):
        self.name = filename.split("/")[-1][:-4]
        self.name = "scRNA"
        self.sequences, self.sequences_by_name = self.get_sequences(filename)
        self.number_sequences = len(self.sequences)

    # returns sequences as list - for iterations - and as dict - for access by name
    @staticmethod
    def get_sequences(filename):
        # loading the file
        raw_data = open(filename, "r").readlines()[1:]
        sequences = []
        sequences_by_name = {}
        for line in raw_data:
            # delimiter in the file is ",", returns list of properties as strings
            line = line.split(",")
            sequences.append(SingleCellRNASequence(line[0], float(line[1]), float(line[2]), float(line[3])))
            sequences_by_name[line[0]] = SingleCellRNASequence(line[0], float(line[1]), float(line[2]), float(line[3]))
        return sequences, sequences_by_name
