#!/usr/bin/env python

import argparse
import re
import cairo

# GLOBAL VARIABLES
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

DRAW_HEIGHT = 5

# ARGUMENT/FILE PARSERS
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
    # parse file
    with open(infile, "r") as rf:
        for line in rf:
            line = line.strip("\n")
            if line.startswith(">"):
                if first_line:
                    index = line
                    first_line = False
                else:
                    record = Gene(index, read)
                    record_list.append(record)
                    index = line
                    read = ''
            else:
                read += line
    record = Gene(index, read)
    record_list.append(record)
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

def gene_splitter(read: str) -> list:
    '''Given a read, splits into a list of features.'''
    # look into re.split for capitals split
    features = []
    working_intron = ''
    working_exon = ''
    caps = read[0].isupper()
    for i, char in enumerate(read):
        # if the character has become upper when the previous was lower, add working_intron to features
        if char.isupper() and caps == False:
            features.append(Feature(working_intron, False))
            working_intron = ''
            caps = True
        # if the character has become lower when the previous was upper, add working_exon to features
        elif char.islower() and caps == True:
            features.append(Feature(working_exon, True))
            working_exon = ''
            caps = False
        # Exon
        if caps:
            working_exon += char
        # Intron
        else:
            working_intron += char
    # add last exon/intron
    if working_intron == '':
        features.append(Feature(working_exon, True))
    else:
        features.append(Feature(working_intron, False))
    return features


# CLASSES
class Gene:
    '''A Gene object. A gene has an index, Features (introns and exons), and a length of the total gene.'''
    def __init__(self, index, seq) -> None:
        '''Pass an index (fasta header or name) and a sequence that is split into a list of Feature objects based on capitalization of the sequence.'''
        self.index = index
        self.length = int(len(seq))
        self.features = gene_splitter(seq)
    
    def __str__(self) -> str:
        '''Print magic method. Returns just the index.'''
        return f"{self.index}"
    
    def __len__(self) -> int:
        '''Length magic method. Returns length of the entire gene, introns and exons.'''
        return self.length

    def count(self) -> int:
        '''Returns integer count of features (introns or exons).'''
        return len(self.features)

    def draw_features(self, surface: cairo.Surface, x: int, y: int):
        '''Given a top-left point, draws a rectangle and line.'''
        context = cairo.Context(surface)
        context.set_line_width(1)
        curr_x = x
        curr_y = y
        # Draw features
        for feature in self.features:
            if feature.is_exon():
                # x,y,width,height 
                context.rectangle(curr_x, curr_y, len(feature), y + DRAW_HEIGHT)
                context.fill()
                curr_x += len(feature)
            else:
                # x,y -> x,y
                context.move_to(curr_x, curr_y + DRAW_HEIGHT)
                context.line_to(curr_x + len(feature), curr_y + DRAW_HEIGHT)
                context.stroke()
                curr_x += len(feature)





class Feature:
    '''A Feature object. A feature is either an exon (exon == True) or an intron (exon == False).'''
    def __init__(self, seq, exon) -> None:
        self.seq = seq
        self.revcomp = revcomp(seq)
        self.length = int(len(seq))
        self.exon = exon
    
    def __str__(self) -> str:
        '''Print magic method.'''
        return self.seq

    def __len__(self) -> int:
        '''Length magic method.'''
        return self.length

    def is_exon(self) -> bool:
        '''Returns True if feature is an exon.'''
        return self.exon


class Motif:
    '''A Motif object.'''
    def __init__(self, motif) -> None:
        '''From passing the motif string, finds attributes length and regex.'''
        self.motif = motif
        self.length = len(motif)
        self.regex = motif_regex(motif)

    def __str__(self) -> str:
        '''Print magic method'''
        return f"{self.motif} with pattern {self.regex}"

    def __len__(self) -> int:
        '''Length magic method.'''
        return self.length


args = get_args()
records = fasta_parser(args.file)
motifs = motif_parser(args.motifs)

testgene = records[0]
print(testgene)

# draw
with cairo.PDFSurface("motif-mark.pdf", 1010, 100) as surface:
    # Make a rectangle for a gene
    print("TESTING DRAW:")
    testgene.draw_features(surface, 5, 5)


    




