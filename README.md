# Bio Linear Aligner
Linear Alignment using BioPython

## How to use

setup venv and install dependencies
```sh
python -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

show all options
```sh
python aligner.py -h
```

```
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
  ```


run local alignment on sample `sample_sequences/two_sequences.fasta` sequences
```sh
python aligner.py local sample_sequences/two_sequences.fasta
```

output
```
Maximum score: 11.0
MKLKKTIGAMALATLFATMGASAVEKTISVTASVDPTVDLLQSDGSALPNSVALTYSPAV                   
                                              ||-|-|-|--|||                    
              MLKIKYLLIGLSLSAMSSYSLAAAGPTLTKELAL-N-V-L--SPAALDATWAPQDNLTLSNTGVS
Score: 11.0
```

you can also provie 2 input files
```sh
python aligner.py global sample_sequences/dna-1.txt sample_sequences/dna-2.txt
```

```output
Maximum score: 5.0
tcgattata----gag
---||-||.----|||
---at-atcccccgag
Score: 5.0

tcgattata----gag
---|-|||.----|||
---a-tatcccccgag
Score: 5.0

tcgattat-a---gag
---||-||-.---|||
---at-atcccccgag
Score: 5.0

tcgattat-a---gag
---|-|||-.---|||
---a-tatcccccgag
Score: 5.0

tcgattat--a--gag
---||-||--.--|||
---at-atcccccgag
Score: 5.0

tcgattat--a--gag
---|-|||--.--|||
---a-tatcccccgag
Score: 5.0

tcgattat---a-gag
---||-||---.-|||
---at-atcccccgag
Score: 5.0

tcgattat---a-gag
---|-|||---.-|||
---a-tatcccccgag
Score: 5.0

tcgattat----agag
---||-||----.|||
---at-atcccccgag
Score: 5.0

tcgattat----agag
---|-|||----.|||
---a-tatcccccgag
Score: 5.0
```

you can customize the scoring and penality with 
```
--matchscore  match score [default: 2]
--mismatchscore  mismatch score [default: -1]
--opengapscore  open gap score [default: -1]
--extendgapscore  extend gap score [default: -1]
```

or you can use one of predefined substitution matrices
["benner6", "benner22", "benner74", "blosum100",
  "blosum30", "blosum35", "blosum40", "blosum45",
  "blosum50", "blosum55", "blosum60", "blosum62",
  "blosum65", "blosum70", "blosum75", "blosum80",
  "blosum85", "blosum90", "blosum95", "feng",
  "fitch", "genetic", "gonnet", "grant",
  "ident", "johnson", "levin", "mclach",
  "miyata", "nwsgappep", "pam120", "pam180",
  "pam250", "pam30", "pam300", "pam60",
  "pam90", "rao", "risler", "structure",
]
```
python aligner.py local sample_sequences/two_sequences.fasta --matrix blosum62
```


default output file is `output.txt` but you can specifiy it with -o OUT_FILE, for example
```sh
python aligner.py local sample_sequences/two_sequences.fasta -o alignments.txt
```

preview
```sh
cat alignments.txt
```