"""
Usage:
  aligner.py (local|global) INPUT ...
  [-o OUT_FILE] 
  [--matchscore MATCH_SCORE] 
  [--mismatchscore MISMATCH_SCORE] 
  [--opengapscore OPEN_GAP_SCORE] 
  [--extendgapscore EXTEND_GAP_SCORE] 
  [--informat FORMAT] 
  [--matrix MATRIX]
  aligner.py (-h | --help)

Examples:
  aligner.py local sample_sequences/two_sequences.fasta
  aligner.py global sample_sequences/dna-1.txt sample_sequences/dna-2.txt
  aligner.py local sample_sequences/two_sequences.fasta -o alignments.txt
  aligner.py local sample_sequences/two_sequences.fasta --matrix blosum62

Options:
  -h --help  Show this screen.
  -o OUT_FILE  output file [default: ./output.txt]
  --matchscore MATCH_SCORE  match score [default: 2]
  --mismatchscore MISMATCH_SCORE  mismatch score [default: -1]
  --opengapscore OPEN_GAP_SCORE  open gap score [default: -1]
  --extendgapscore EXTEND_GAP_SCORE  extend gap score [default: -1]
  --informat FORMAT  input format [default: fasta]
  --matrix MATRIX  substitution matrice
"""

from docopt import docopt
from Bio.SeqIO import parse, read
from Bio.Data import IUPACData
from Bio import SeqIO, AlignIO, pairwise2, Align
from Bio.Align import substitution_matrices
from Bio.Alphabet import generic_dna, generic_protein, generic_rna, IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, transcribe, back_transcribe
from Bio.SeqUtils import GC
from Bio.Data import CodonTable
from Bio.SubsMat import MatrixInfo

if __name__ == '__main__':
    arguments = docopt(__doc__)
    input_files = arguments['INPUT']
    aligner = Align.PairwiseAligner()
    if arguments['local']:
        aligner.mode = 'local'
    else:
        aligner.mode = 'global'

    if arguments['--matrix'] in MatrixInfo.available_matrices:
        matrice = arguments['--matrix']
        loaded_matrix = substitution_matrices.load(matrice.upper())
        aligner.substitution_matrix =  loaded_matrix
    elif arguments['--matrix']:
        matrix_file = arguments['--matrix']
        with open(matrix_file, 'r') as matrix_file:
            matrice_dict = substitution_matrices.read(matrix_file)
            aligner.substitution_matrix = matrice_dict
    else:
        aligner.match_score = int(arguments['--matchscore'])
        aligner.mismatch_score = int(arguments['--mismatchscore'])
        aligner.open_gap_score = int(arguments['--opengapscore'])
        aligner.extend_gap_score = int(arguments['--extendgapscore'])

    out_file = arguments['-o']
    input_file_format = arguments['--informat']

    if len(input_files) == 2:
        sequences = [read(input_file, input_file_format)
                     for input_file in input_files]
        first = sequences[0]
        second = sequences[1]
    elif len(input_files) == 1:
        sequences = parse(input_files[0], input_file_format)
        first = next(sequences)
        second = next(sequences)
    else:
        raise(RuntimeError('Please provide exactly 1 or 2 sequences in each file'))

    score = aligner.score(first.seq, second.seq)
    print(f'Maximum score: {score}')
    alignments = aligner.align(first.seq, second.seq)
    with open(out_file, 'w') as out_file:
        for align in alignments:
            print(align.format(), end="")
            print(f'Score: {align.score}', end="\n\n")
            print(align.format(), end="", file=out_file)
            print(f'Score: {align.score}', end="\n\n", file=out_file)

