#!/usr/bin/env python

from pathlib import Path
import os
import argparse

import filter_msa
import check_complex

dir_temp = Path("temp")
dir_temp.mkdir(parents=True, exist_ok=True)

dir_result = Path("result")
dir_result.mkdir(parents=True, exist_ok=True)

def string_to_bool(value):
    if value.lower() in {"false", "f", "no", "n", "0"}:
        return False
    elif value.lower() in {"true", "t", "yes", "y", "1"}:
        return True
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


parser = argparse.ArgumentParser()
parser.add_argument("-m", required=True, help="msa file to give (either .fasta or .a3m)", type=str)
parser.add_argument("-s", required=True, help="structure file (either .pdb or .cif)", type=str)
parser.add_argument("-i", required=True, help="maximum percentage identity between sequences to keep", type=int)
parser.add_argument("-f", required=True, help="final number of sequences to keep", type=int)
parser.add_argument("--split", type=string_to_bool, required=True, help="split fasta before filtering sequences")


args = parser.parse_args()

input_msa = args.m
# input_msa = "complex/2_HBA_2_HBB_HUMAN.a3m"

structure = args.s
A3M_CONVERTER = "/data/work/I2BC/hugo.pointier/msa_tools/script/a3m_to_fasta.sh"

# TODO: delete
# A3M_CONVERTER = "./a3m_to_fasta.sh"


def check_msa(input):
    """check if input is fasta or a3m and convert it if it is the case"""
    global dir_temp
    fasta_msa = ""
    valid_extension = [".fasta", ".a3m"]
    if Path(input).is_file() and Path(input).suffix in valid_extension:
        if Path(input).suffix == ".fasta":
            fasta_msa = input_msa
        else:
            fasta_msa = str(dir_temp) + "/" + str(Path(input).stem) + ".fasta"
            os.system(A3M_CONVERTER + f" -i {input_msa} -o {fasta_msa}")

    else:
        print("Wrong filetype")
    return fasta_msa


fasta_msa = check_msa(input_msa)
obj = check_complex.r4s_multi(
    msa_input=fasta_msa,
    output="output",
    structure_file=structure,
    dir_temp=str(dir_temp),
    maximum_percentage_identity=args.i,
    maximum_number_sequences=args.f,
    split_identity=args.split,
)
obj.run()
