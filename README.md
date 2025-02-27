# motif-mark

A tool for locating and visualizing splicing motifs on sequences. See usage for file requirements. Outputs an annotated PNG output based on input FASTA name.

## Usage:

`motif-mark-oop.py [-h] -f FASTA -m MOTIF`

FASTA file with max 10 sequences and ≤1000 bases per sequence.
Motifs file with one motif per line in a text file and <≤10 bases per motif. Bases must comply with [IUPAC naming conventions](https://genome.ucsc.edu/goldenPath/help/iupac.html).

## Requirements
[Pycairo v.1.27.0](https://github.com/pygobject/pycairo)

## Example Output
Example here
