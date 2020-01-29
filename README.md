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


default output file is `output.txt` but you can specifiy it with -o OUT_FILE, for example
```sh
python aligner.py local sample_sequences/two_sequences.fasta -o alignments.txt
```

preview
```sh
cat alignments.txt
```