"""
Usage:
  aligner.py (local|global) INPUT ...
  [-o OUT_FILE] 
  [--matchscore match_score] 
  [--mismatchscore mismatch_score] 
  [--opengapscore open_gap_score] 
  [--extendgapscore extend_gap_score] 
  [--informat FORMAT] 
  aligner.py (-h | --help)

Examples:
  aligner.py local sample_sequences/two_sequences.fasta
  aligner.py global sample_sequences/dna-1.txt sample_sequences/dna-2.txt
  aligner.py local sample_sequences/two_sequences.fasta -o alignments.txt

Options:
  -h --help  Show this screen.
  -o OUT_FILE  Output file [default: ./output.txt]
  --matchscore  match score [default: 2]
  --mismatchscore  mismatch score [default: -1]
  --opengapscore  open gap score [default: -1]
  --extendgapscore  extend gap score [default: -1]
  --format  input format [default: fasta]
"""

from docopt import docopt
from Bio.SeqIO import parse, read
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
    print(arguments)

    input_files = arguments['INPUT']
    print(input_files)

    aligner = Align.PairwiseAligner()
    if arguments['local']:
        aligner.mode = 'local'
    else:
        aligner.mode = 'global'

    aligner.match_score = arguments['match_score']
    aligner.mismatch_score = arguments['mismatch_score']
    aligner.open_gap_score = arguments['open_gap_score']
    aligner.extend_gap_score = arguments['extend_gap_score']

    out_file = arguments['-o']
    input_file_format = arguments['FORMAT']

    if len(input_files) == 2:
        sequences = [read(input_file, input_file_format)
                     for input_file in input_files]
        first = sequences[0]
        second = sequences[1]
        print(first.seq)
        print(second.seq)
    elif len(input_files) == 1:
        sequences = parse(input_files[0], input_file_format)
        first = next(sequences)
        second = next(sequences)
        print(first.seq)
        print(second.seq)
    else:
        raise(RuntimeError('Please provide exactly 1 or 2 sequences in each file'))

    print(aligner)
    score = aligner.score(first.seq, second.seq)
    print(f'Maximum score: {score}')
    alignments = aligner.align(first.seq, second.seq)
    with open(out_file, 'w') as out_file:
        for align in alignments:
            print(align.format(), end="")
            print(f'Score: {align.score}', end="\n\n")
            print(align.format(), end="", file=out_file)
            print(f'Score: {align.score}', end="\n\n", file=out_file)
