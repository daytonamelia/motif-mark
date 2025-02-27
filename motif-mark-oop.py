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

distinct_colors = {
    "Maroon": [128,0,0],
    "Brown": [170,110,40],
    "Teal": [0,128,128],
    "Navy": [0,0,128],
    "Red": [230,25,75],
    "Orange": [245,130,48],
    "Yellow": [255,255,25],
    "Green": [60,180,75],
    "Cyan": [70,240,240],
    "Blue": [0,130,200],
    "Magenta": [240,50,230],
    "Grey": [128,128,128],
    "Pink": [250,190,212],
    "Mint": [170,255,195],
    "Lavender": [220,190,255]}

MARGIN = 10
SPACING = 10
DRAW_HEIGHT = 10
FONT_SIZE = 5
FONT_FACE = "Arial"

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
    colors = ["Blue", "Orange", "Maroon", "Lavender", "Navy"]
    with open(infile, "r") as rf:
        for i, line in enumerate(rf):
            line = line.strip("\n")
            motif = Motif(line, colors[i])
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
        self.seq = seq
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

    def draw_features(self, surface: cairo.Surface, surface_x: int, surface_y: int):
        '''Given a top-left point, writes the index, then draws a rectangle and line.'''
        # Context and variables
        context = cairo.Context(surface)
        context.set_line_width(1)
        curr_x = surface_x
        curr_y = surface_y
        # Write the index
        context.set_font_size(FONT_SIZE)
        context.select_font_face(FONT_FACE)
        context.move_to(curr_x, curr_y) # x,y
        context.show_text(self.index)
        context.stroke()
        curr_y += FONT_SIZE
        # Draw the features
        for feature in self.features:
            if feature.is_exon(): # draw exons 
                context.rectangle(curr_x, curr_y, len(feature), DRAW_HEIGHT * 2) # x,y,width,height
                context.fill()
                curr_x += len(feature)
            else:  # draw introns
                context.move_to(curr_x, curr_y + DRAW_HEIGHT) # x,y
                context.line_to(curr_x + len(feature), curr_y + DRAW_HEIGHT) #x,y
                context.stroke()
                curr_x += len(feature)


class Feature:
    '''A Feature object. A feature is either an exon (exon == True) or an intron (exon == False) that might include motifs.'''
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
    '''A Motif object, with a motif, a regex pattern, and a color for drawing.'''
    def __init__(self, motif, color) -> None:
        '''From passing the motif string, finds attributes length and regex.'''
        self.motif = motif
        self.length = len(motif)
        self.regex = motif_regex(motif)
        self.color = distinct_colors[color]

    def __str__(self) -> str:
        '''Print magic method'''
        return self.motif

    def __len__(self) -> int:
        '''Length magic method.'''
        return self.length

        self.color = color


# Set up and parsing
args = get_args()
records = fasta_parser(args.file)
motifs = motif_parser(args.motifs)

# Draw features
with cairo.PDFSurface("motif-mark.pdf", 1010, 40 * len(records) + 35) as surface:
    # Context variables and surface coordinates setup
    surface_x = MARGIN
    surface_y = MARGIN
    context = cairo.Context(surface)
    # Make legend colors and text
    context.set_font_size(FONT_SIZE)
    context.select_font_face(FONT_FACE)
    for motif in motifs:
        # draw color box
        context.set_line_width(1)
        context.set_source_rgb(motif.color[0]/255, motif.color[1]/255, motif.color[2]/255)
        context.rectangle(surface_x, surface_y, FONT_SIZE, FONT_SIZE)
        context.fill()
        # draw text
        context.set_line_width(1)
        context.set_source_rgb(0,0,0)
        context.move_to(surface_x + FONT_SIZE * 2, surface_y + FONT_SIZE/1.1) # x,y
        context.show_text(motif.motif)
        context.stroke()
        surface_y += FONT_SIZE + SPACING/2
    # Make legend box
    context.rectangle(MARGIN/1.5, MARGIN/2, 100, surface_y - FONT_SIZE)
    context.stroke()
    # Add spacing
    surface_y += SPACING
    # Draw genes
    for feature in records:
        feature.draw_features(surface, surface_x, surface_y)
        surface_y += FONT_SIZE + DRAW_HEIGHT + MARGIN + SPACING


    




