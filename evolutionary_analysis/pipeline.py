#!/usr/bin/env python

from pathlib import Path
import os
import argparse

import filter_msa
import check_complex

dir_temp = Path("temp")
dir_temp.mkdir(parents=True, exist_ok=True)


def string_to_bool(value):
    if value.lower() in {"false", "f", "no", "n", "0"}:
        return False
    elif value.lower() in {"true", "t", "yes", "y", "1"}:
        return True
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


parser = argparse.ArgumentParser()
parser.add_argument(
    "-m",
    required=True,
    nargs="+",
    help="msa file to give (either .fasta or .a3m)",
    type=str,
)
parser.add_argument(
    "-s", required=True, help="structure file (either .pdb or .cif)", type=str
)
parser.add_argument(
    "-i",
    required=True,
    help="maximum percentage identity between sequences to keep",
    type=int,
)
parser.add_argument(
    "-f", required=True, help="final number of sequences to keep", type=int
)
parser.add_argument("-o", required=True, help="output directory", type=str)
parser.add_argument(
    "--split",
    required=True,
    help="split fasta before filtering sequences",
    type=string_to_bool,
)
parser.add_argument(
    "--chains", required=False, nargs="+", help="array of chain wanted", type=str
)
parser.add_argument(
    "--multiple",
    required=True,
    help="if there is multiple msa file or not",
    type=string_to_bool,
)
parser.add_argument(
    "--field",
    required=False,
    help="occupancy / bfactor => in which column should the conservation be added",
    type=str,
    default="occupancy",
)

args = parser.parse_args()

input_msa = args.m

structure = args.s
A3M_CONVERTER = os.path.join(os.path.dirname(os.path.abspath(__file__)),"a3m_to_fasta.sh")

# TODO: delete
# A3M_CONVERTER = "./a3m_to_fasta.sh"


def check_msa(input, output_directory):
    """check if input is fasta or a3m and convert it if it is the case"""
    if not isinstance(input,list):
        input = [input]
    fasta_msas = []
    for myinput in input:
        fasta_msa = ""
        valid_extension = [".fasta", ".fa", ".a3m"]
        if Path(myinput).is_file() and Path(myinput).suffix in valid_extension:
            if Path(myinput).suffix == ".fasta" or Path(myinput).suffix == ".fa":
                fasta_msa = myinput
            else:
                fasta_msa = os.path.join(output_directory, str(Path(myinput).stem) + ".fasta")
                os.system(A3M_CONVERTER + f" -i {myinput} -o {fasta_msa}")
            fasta_msas.append(fasta_msa)
        else:
            print("Wrong filetype")
    return fasta_msas[0] if (not isinstance(input,list)) or (isinstance(input,list) and len(input)==1) else fasta_msas


dir_result = Path(args.o)
dir_result.mkdir(parents=True, exist_ok=True)

if not args.multiple:
    fasta_msa = check_msa(input_msa, dir_temp)
else:
    fasta_msa = []
    for file in input_msa:
        fasta_msa.append(check_msa(file, dir_temp))

obj = check_complex.r4s_multi(
    msa_input=fasta_msa,
    output_directory=args.o,
    array_chains=None if not args.chains else [x.upper() for x in args.chains],
    structure_file=structure,
    temporary_directory=str(dir_temp),
    maximum_percentage_identity=args.i,
    maximum_number_sequences=args.f,
    split_identity=args.split,
    multiple_msa=args.multiple,
    field=args.field,
)
obj.run()
