"""
Usage:
  aligner.py (local|global) 
    [--matchscore match_score]
    [--mismatchscore mismatch_score]
    [--opengapscore open_gap_score] 
    [--extendgapscore extend_gap_score] 
    [FILE] [--format FORMAT] 
    [-o OUT_FILE]
  aligner.py (-h | --help)

Options:
  -h --help         Show this screen.
  --matchscore      Default 2,
  --mismatchscore   Default -1,
  --opengapscore    Default -1,
  --extendgapscore  Default -1,
  --format          Default fasta,
"""

from docopt import docopt
from Bio.SeqIO import parse
from Bio.Data import IUPACData
from Bio import SeqIO, AlignIO, pairwise2, Align 
from Bio.Alphabet import generic_dna, generic_protein, generic_rna, IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, transcribe, back_transcribe
from Bio.SeqUtils import GC
from Bio.Data import CodonTable 
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

if __name__ == '__main__':
    arguments = docopt(__doc__)

    if not arguments['FILE']:
        input_file = "input.fasta"
    else:
        input_file = arguments['FILE']

    if not arguments['OUT_FILE']:
        out_file = "alignments.txt"
    else:
        out_file = arguments['OUT_FILE']

    if not arguments['FORMAT']:
        input_file_format = "fasta"
    else:
        input_file_format = arguments['FORMAT']

    aligner = Align.PairwiseAligner()
    if arguments['local']:
        aligner.mode = 'local'
    else:
        aligner.mode = 'global'
    
    if arguments['match_score']:
        aligner.match_score = arguments['match_score']
    else:
        aligner.match_score = 2
    
    if arguments['mismatch_score']:
        aligner.mismatch_score = arguments['mismatch_score']
    else:
        aligner.mismatch_score = -1

    if arguments['open_gap_score']:
        aligner.open_gap_score = arguments['open_gap_score']
    else:
        aligner.open_gap_score = -1
    
    if arguments['extend_gap_score']:
        aligner.extend_gap_score = arguments['extend_gap_score']
    else:
        aligner.extend_gap_score = -1


    sequences = parse(input_file, input_file_format)
    first = next(sequences)
    print(first.seq)
    second = next(sequences)
    print(second.seq)

    print(aligner)
    score = aligner.score(first.seq, second.seq)
    print(f'Maximum score: {score}')
    alignments = aligner.align(first.seq, second.seq)
    with open(out_file, 'w') as out_file:
        for align in alignments:
            print(align.format(), end="", file=out_file)
            print(f'Score: {align.score}', end="\n\n", file=out_file)
