# Bio Linear Aligner
Linear Alignment using BioPython

## How to use

setup venv and install dependencies
```sh
python -m venv venv
. venv/bin/activate
pip insall -r requirements.txt
```

show all options
```sh
python aligner.py -h
```

run local alignment on sample `input.fasta` sequences
```sh
python aligner.py local
```

preview
```sh
cat alignments.txt
```