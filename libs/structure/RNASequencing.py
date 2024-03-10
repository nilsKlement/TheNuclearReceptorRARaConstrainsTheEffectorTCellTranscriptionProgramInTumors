import os, sys
sys.path.insert(1, os.getcwd())
sys.path.insert(1, os.getcwd()+"/libs/structure")
from RNASequence import RNASequence, RNASequenceGeneral


class RNASequencingGeneral:

    # a class based on the general RNA sequence for a RNA sequencing

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
            sequences.append(RNASequenceGeneral(line[0], line[1:]))
            sequences_by_name[line[0]] = RNASequenceGeneral(line[0], line[1:])
        return sequences, sequences_by_name


class RNASequencing:

    # a class based on the specific RNA sequence for a RNA sequencing

    def __init__(self, filename):
        self.name = filename.split("/")[-1][:-4]
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
            # removing incomplete or false data
            if line[-1] != "\n" and line[-2] != '#DIV/0!' and line[-1] != '#DIV/0!\n' and line[-2] != '#NUM!' and line[-1] != '#NUM!':
                # appending specific RNA sequence
                sequences.append(RNASequence(line[-3], float(line[-2]), float(line[-1])))
                sequences_by_name[line[-3]] = RNASequence(line[-3], float(line[-2]), float(line[-1]))
        return sequences, sequences_by_name
