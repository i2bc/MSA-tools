#!/usr/bin/env python

import os

HHFILTER = "~/programmation/stage/script_msa_tools/hhsuite/bin/hhfilter"


def run_hhfilter(input_file, output_file, max_percentage_identity=80):
    os.system(f"{HHFILTER} -i {input_file} -o {output_file} -id {max_percentage_identity}")


def select_remaining(input_file, output_file, maximum_sequences):
    with open(input_file, "r") as f:
        with open(output_file, "w") as o:
            lines = f.readlines()
            counter = 0
            for line in lines:
                if line.startswith(">"):
                    counter += 1
                if counter > maximum_sequences:
                    break
                o.write(line)
