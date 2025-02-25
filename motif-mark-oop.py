#!/usr/bin/env python

import argparse
import re
# import cairo

# Global variables
ambig_nts = {
    "R": "(A|G)",
    "K": "(G|T)",
    "S": "(G|C)",
    "Y": "(C|T)",
    "M": "(A|C)",
    "W": "(A|T)",
    "B": "[^A]",
    "H": "[^G]",
    "D": "[^C]",
    "V": "[^T]",
    "N": "(A|T|C|G)"}


# Argument/file parsers
def get_args():
    parser = argparse.ArgumentParser(description="Visualize motifs on sequences.")
    parser.add_argument("-f", "--file", help="Designates absolute file path to fasta file", type=str, required=True)
    parser.add_argument("-m", "--motifs", help="Designates absolute file path to motif file (one motif per line)", type=str, required=True)
    return parser.parse_args()

def fasta_parser(infile:str) -> list:
    '''Given a fasta file with reads that might span multiple lines, return a list of Record objects.'''
    record_list = []
    index = ''
    read = ''
    first_line = True
    with open(infile, "r") as rf:
        for line in rf:
            line = line.strip("\n")
            if line.startswith(">"):
                if first_line:
                    index = line
                    first_line = False
                else:
                    record = Record(index, read)
                    record_list.append(record)
                    index = line
            else:
                read += line
    return record_list

def motif_parser(infile:str) -> list:
    '''Given a motif file of one motif per line, return a list of Motif objects.'''
    motif_list = []
    with open(infile, "r") as rf:
        for line in rf:
            line = line.strip("\n")
            motif = Motif(line)
            motif_list.append(motif)
    return motif_list


# SIMPLE FUNCTIONS
def revcomp(DNA:str) -> str:
    '''Returns the reverse complement of a DNA sequence.'''
    DNAtable = str.maketrans("ATCG", "TAGC")
    return DNA[::-1].translate(DNAtable)


# COMPLEX FUNCTIONS
def motif_regex(motif:str) -> re.Pattern:
    '''Given a motif string, transform into a regex object. Requires the python re module.'''
    regex_str = ''
    motif = motif.strip()
    for char in motif:
        char = char.upper()
        if char in ["A", "T", "C", "G"]:
            regex_str += char
        elif char == "U":
            regex_str += "T"
        else:
            regex_str += ambig_nts[char]
    return re.compile(regex_str)

def gene_splitter(read: str) -> list, list:
    '''Given a read, splits into a list of introns and exons.'''
    # look into re.split for capitals split
    pass

# CLASSES
class Record:
    def __init__(self, index, read):
        self.index = index
        self.read = read
        self.revcomp = revcomp(read)
        self.length = len(read)
        
        self.introns = []
        self.exons = []
    
    def __str__(self):
        '''Print magic method.'''
        return f"{self.index}\n{self.read}"
    
    def __len__(self):
        return self.length

class Motif:
    '''TODO'''
    def __init__(self, motif):
        '''From passing the motif string, finds attributes length and regex.'''
        self.motif = motif
        self.length = len(motif)
        self.regex = motif_regex(motif)

    def __str__(self):
        '''Print magic method'''
        return f"{self.motif} with pattern {self.regex}"

    def __len__(self):
        '''Length magic method.'''
        return self.length


args = get_args()
records = fasta_parser(args.file)
motifs = motif_parser(args.motifs)






