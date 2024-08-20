import os
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, help="input file, either a3m or fasta")

args = parser.parse_args()

input = args.i

extension = Path(input).suffix
filename_without_extension = Path(input).stem

A3M_TO_FASTA = "./a3m_to_fasta.sh"
FASTA_TO_A3M = "./fasta_to_a3m.sh"

# the validator in django only accept .a3m, .fa or .fasta
if extension == ".a3m":
    os.system(A3M_TO_FASTA + f" -i {input} -o {filename_without_extension}.fasta")
else:
    os.system(FASTA_TO_A3M + f" -i {input} -o {filename_without_extension}.a3m")
