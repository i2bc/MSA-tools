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

dir_zip = Path("zip")
dir_zip.mkdir(parents=True, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument("-m", required=True)
parser.add_argument("-s", required=True)
parser.add_argument("-i", required=True)
parser.add_argument("-f", required=True)

args = parser.parse_args()

input_msa = args.m
# input_msa = "complex/2_HBA_2_HBB_HUMAN.a3m"

structure = args.s
A3M_CONVERTER = "/data/work/I2BC/hugo.pointier/msa_tools/script/a3m_to_fasta.sh"


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


def filter_fasta(
    fasta_msa,
    hh_filter_msa,
    filtered_msa,
    max_percentage_identity=80,
    maximum_sequences=100,
):
    filter_msa.run_hhfilter(fasta_msa, hh_filter_msa, max_percentage_identity)
    filter_msa.select_remaining(hh_filter_msa, filtered_msa, maximum_sequences)


hh_filter = os.path.join(str(dir_temp), "hh_filter.fasta")
filtered_msa = os.path.join(str(dir_temp), "filtered.fasta")
output = os.path.join(dir_result, "output")
fasta_msa = check_msa(input_msa)
filter_fasta(fasta_msa, hh_filter, filtered_msa, int(args.i), int(args.f))

obj = check_complex.r4s_multi(
    msa_input=filtered_msa,
    output=output,
    structure_file=structure,
    dir_temp=str(dir_temp),
)
obj.run()
